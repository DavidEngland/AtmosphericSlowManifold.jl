function _cal_feature_key(f::StateVariable)
    return f.name
end

function _cal_feature_key(f::DiagnosticVariable)
    return f.name
end

function _cal_finite_diff_col(col::Vector{Float64}, z::Vector{Float64}, order::Int)
    order >= 1 || return copy(col)
    d = copy(col)
    for _ in 1:order
        out = zeros(length(col))
        out[1] = (d[2] - d[1]) / (z[2] - z[1])
        for i in 2:(length(col) - 1)
            out[i] = (d[i + 1] - d[i - 1]) / (z[i + 1] - z[i - 1])
        end
        out[end] = (d[end] - d[end - 1]) / (z[end] - z[end - 1])
        d = out
    end
    return d
end

function _cal_feature_column(obs::ObservationTable, feature::AbstractBasisFeature)
    if feature isa StateVariable || feature isa DiagnosticVariable
        key = _cal_feature_key(feature)
        haskey(obs.columns, key) || throw(KeyError(key))
        return Float64.(obs.columns[key])
    elseif feature isa SpatialDerivative
        haskey(obs.columns, :z) || throw(ArgumentError("ObservationTable must include :z for derivative features."))
        key = feature.variable
        haskey(obs.columns, key) || throw(KeyError(key))
        return _cal_finite_diff_col(Float64.(obs.columns[key]), Float64.(obs.columns[:z]), feature.order)
    else
        throw(ArgumentError("Unsupported feature type $(typeof(feature))."))
    end
end

function _cal_term_column(obs::ObservationTable, term::OperatorTerm)
    n = length(first(values(obs.columns)))
    isempty(term.basis) && return ones(n)

    col = ones(Float64, n)
    for b in term.basis
        base = _cal_feature_column(obs, b.feature)
        col .*= base .^ b.power
    end
    return col
end

function _cal_target_symbol(model::DiscoveredModel, obs::ObservationTable)
    if haskey(obs.columns, model.target_variable)
        return model.target_variable
    elseif haskey(obs.columns, :u)
        return :u
    end
    return first(keys(obs.columns))
end

function _cal_design_matrix(obs::ObservationTable, model::DiscoveredModel)
    n_terms = length(model.terms)
    n_terms == 0 && return zeros(length(first(values(obs.columns))), 0)
    cols = [_cal_term_column(obs, t) for t in model.terms]
    return hcat(cols...)
end

function _cal_fit_ridge(X::Matrix{Float64}, y::Vector{Float64}, λ::Float64)
    p = size(X, 2)
    p == 0 && return Float64[]
    A = X' * X + λ * I(p)
    b = X' * y
    return A \ b
end

function _cal_update_model(model::DiscoveredModel{Float64}, θ::Vector{Float64}, residual_norm::Float64)
    length(θ) == length(model.terms) || throw(ArgumentError("Coefficient length mismatch during model update."))
    new_terms = OperatorTerm{Float64}[]
    for (i, t) in enumerate(model.terms)
        push!(new_terms, OperatorTerm{Float64}(θ[i], t.basis))
    end
    sparsity = count(abs.(θ) .> 1e-10)
    return DiscoveredModel{Float64}(model.target_variable, new_terms, residual_norm, sparsity)
end

function dispatch_calibrate(
    model::DiscoveredModel{Float64},
    obs::ObservationTable,
    alg::MaximumLikelihood,
)
    ysym = _cal_target_symbol(model, obs)
    y = Float64.(obs.columns[ysym])
    X = _cal_design_matrix(obs, model)

    λ = max(alg.tol, 1e-12)
    mle = _cal_fit_ridge(X, y, λ)
    yhat = isempty(mle) ? zeros(length(y)) : X * mle
    residual = yhat .- y
    rss = dot(residual, residual)
    ybar = sum(y) / max(length(y), 1)
    tss = sum((y .- ybar) .^ 2)
    r2 = tss > 0 ? 1.0 - rss / tss : 1.0

    calibrated_model = _cal_update_model(model, mle, norm(residual))

    params = Dict{Symbol, Vector{Float64}}(:mle => mle)
    diagnostics = Dict{Symbol, Any}(
        :status => :ok,
        :backend => :MaximumLikelihood,
        :optimizer => alg.optimizer,
        :max_iter => alg.max_iter,
        :target_symbol => ysym,
        :residual_rss => rss,
        :r2 => r2,
        :n_obs => length(first(values(obs.columns))),
    )
    return CalibrationResult{MaximumLikelihood}(calibrated_model, params, diagnostics, alg)
end
