# src/Calibration/UncertaintyPropagation.jl
using Statistics

export evaluate_profile_uncertainty

function _up_feature_profile(
    f::StateVariable,
    z_grid::AbstractVector{<:Real},
    feature_profiles::AbstractDict{Symbol, Any},
    scalar_features::AbstractDict{Symbol, <:Real},
)
    if f.name == :z
        return Float64.(z_grid)
    elseif haskey(feature_profiles, f.name)
        v = Float64.(feature_profiles[f.name])
        length(v) == length(z_grid) || throw(ArgumentError("Feature profile $(f.name) length must match z_grid length."))
        return v
    elseif haskey(scalar_features, f.name)
        return fill(Float64(scalar_features[f.name]), length(z_grid))
    end
    throw(ArgumentError("Missing feature profile for state variable $(f.name). Provide it via feature_profiles or scalar_features."))
end

function _up_feature_profile(
    f::DiagnosticVariable,
    z_grid::AbstractVector{<:Real},
    feature_profiles::AbstractDict{Symbol, Any},
    scalar_features::AbstractDict{Symbol, <:Real},
)
    if f.name == :z
        return Float64.(z_grid)
    elseif haskey(feature_profiles, f.name)
        v = Float64.(feature_profiles[f.name])
        length(v) == length(z_grid) || throw(ArgumentError("Feature profile $(f.name) length must match z_grid length."))
        return v
    elseif haskey(scalar_features, f.name)
        return fill(Float64(scalar_features[f.name]), length(z_grid))
    end
    throw(ArgumentError("Missing feature profile for diagnostic variable $(f.name). Provide it via feature_profiles or scalar_features."))
end

function _up_feature_profile(
    f::SpatialDerivative,
    z_grid::AbstractVector{<:Real},
    feature_profiles::AbstractDict{Symbol, Any},
    scalar_features::AbstractDict{Symbol, <:Real},
)
    key = Symbol("d$(f.order)_$(f.variable)")
    if haskey(feature_profiles, key)
        v = Float64.(feature_profiles[key])
        length(v) == length(z_grid) || throw(ArgumentError("Feature profile $(key) length must match z_grid length."))
        return v
    elseif haskey(scalar_features, key)
        return fill(Float64(scalar_features[key]), length(z_grid))
    end
    throw(ArgumentError("Missing feature profile for spatial derivative $(key). Provide it via feature_profiles or scalar_features."))
end

function _up_draw_matrix(model::DiscoveredModel, posterior_draws::AbstractMatrix{<:Real})
    n_terms = length(model.terms)
    n_terms == 0 && return zeros(0, size(posterior_draws, 2))

    if size(posterior_draws, 1) == n_terms
        return Float64.(posterior_draws)
    elseif size(posterior_draws, 2) == n_terms
        return permutedims(Float64.(posterior_draws), (2, 1))
    end

    throw(ArgumentError("posterior_draws must have one dimension equal to n_terms=$(n_terms)."))
end

function _up_draw_matrix(model::DiscoveredModel, posterior_draws::AbstractDict{Symbol, Any})
    n_terms = length(model.terms)
    n_terms == 0 && return zeros(0, 0)

    if haskey(posterior_draws, :draws)
        draws = posterior_draws[:draws]
        draws isa AbstractMatrix || throw(ArgumentError("posterior_draws[:draws] must be a matrix."))
        return _up_draw_matrix(model, draws)
    end

    keys_needed = [Symbol("c_", i) for i in 1:n_terms]
    for k in keys_needed
        haskey(posterior_draws, k) || throw(ArgumentError("Missing posterior draw key $(k)."))
    end

    ndraws = length(posterior_draws[keys_needed[1]])
    ndraws > 0 || throw(ArgumentError("Posterior draw vectors cannot be empty."))

    mat = zeros(n_terms, ndraws)
    for (i, k) in enumerate(keys_needed)
        vals = Float64.(posterior_draws[k])
        length(vals) == ndraws || throw(ArgumentError("All posterior draw vectors must have the same length."))
        mat[i, :] = vals
    end
    return mat
end

function _build_profile_design(
    model::DiscoveredModel,
    z_grid::AbstractVector{<:Real},
    feature_profiles::AbstractDict{Symbol, Any},
    scalar_features::AbstractDict{Symbol, <:Real},
)
    nz = length(z_grid)
    n_terms = length(model.terms)
    n_terms == 0 && return zeros(nz, 0)

    X = ones(Float64, nz, n_terms)

    for (j, term) in enumerate(model.terms)
        if isempty(term.basis)
            X[:, j] .= 1.0
            continue
        end

        col = ones(Float64, nz)
        for b in term.basis
            base = _up_feature_profile(b.feature, z_grid, feature_profiles, scalar_features)
            col .*= base .^ b.power
        end
        X[:, j] = col
    end
    return X
end

"""
    evaluate_profile_uncertainty(model, posterior_draws, z_grid; feature_profiles=Dict(), scalar_features=Dict())

Evaluate sampled eddy diffusivity profiles K_m(z | theta_i) over posterior draws and
return median, 2.5%, and 97.5% quantile ribbons.
"""
function evaluate_profile_uncertainty(
    model::DiscoveredModel,
    posterior_draws::Union{AbstractMatrix{<:Real}, AbstractDict{Symbol, Any}},
    z_grid::AbstractVector{<:Real};
    feature_profiles::AbstractDict{Symbol, Any} = Dict{Symbol, Any}(),
    scalar_features::AbstractDict{Symbol, <:Real} = Dict{Symbol, Float64}(),
)
    length(z_grid) > 0 || throw(ArgumentError("z_grid cannot be empty."))

    draws = _up_draw_matrix(model, posterior_draws)
    X = _build_profile_design(model, z_grid, feature_profiles, scalar_features)

    n_terms = size(X, 2)
    size(draws, 1) == n_terms || throw(ArgumentError("Posterior draw coefficient dimension does not match model term count."))

    if n_terms == 0
        nz = length(z_grid)
        profiles = zeros(nz, max(size(draws, 2), 1))
        zeros_profile = zeros(nz)
        return (
            profiles = profiles,
            median = zeros_profile,
            lower = zeros_profile,
            upper = zeros_profile,
            quantiles = (0.025, 0.5, 0.975),
        )
    end

    profiles = X * draws
    nz, ns = size(profiles)
    ns > 0 || throw(ArgumentError("No posterior samples were provided."))

    med = zeros(nz)
    low = zeros(nz)
    up = zeros(nz)

    for i in 1:nz
        row = vec(@view profiles[i, :])
        low[i] = quantile(row, 0.025)
        med[i] = quantile(row, 0.5)
        up[i] = quantile(row, 0.975)
    end

    return (
        profiles = profiles,
        median = med,
        lower = low,
        upper = up,
        quantiles = (0.025, 0.5, 0.975),
    )
end
