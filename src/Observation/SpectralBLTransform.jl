function _obs_gegenbauerC(n::Int, lambda::Float64, x::Float64)
    n == 0 && return 1.0
    n == 1 && return 2.0 * lambda * x

    c_nm2 = 1.0
    c_nm1 = 2.0 * lambda * x
    for k in 2:n
        c_n = (2.0 * (k + lambda - 1.0) * x * c_nm1 - (k + 2.0 * lambda - 2.0) * c_nm2) / k
        c_nm2 = c_nm1
        c_nm1 = c_n
    end
    return c_nm1
end

function _safe_xcoord(z::Vector{Float64})
    zmin, zmax = extrema(z)
    zmax > zmin || throw(ArgumentError("Observation z-range must be non-degenerate for projection."))
    return @. 2.0 * (z - zmin) / (zmax - zmin) - 1.0
end

function _trapz_inner(x::Vector{Float64}, f::Vector{Float64}, g::Vector{Float64}, lambda::Float64)
    acc = 0.0
    for i in 1:(length(x) - 1)
        x1 = x[i]
        x2 = x[i + 1]
        w1 = max(1e-12, (1.0 - x1^2)^(lambda - 0.5))
        w2 = max(1e-12, (1.0 - x2^2)^(lambda - 0.5))
        y1 = f[i] * g[i] * w1
        y2 = f[i + 1] * g[i + 1] * w2
        acc += 0.5 * (y1 + y2) * (x2 - x1)
    end
    return acc
end

function _project_profile_to_basis(z::Vector{Float64}, y::Vector{Float64}; n_modes::Int, lambda::Float64)
    x = _safe_xcoord(z)
    modes = [Float64[_obs_gegenbauerC(n - 1, lambda, xq) for xq in x] for n in 1:n_modes]

    coeffs = zeros(n_modes)
    for n in 1:n_modes
        numer = _trapz_inner(x, y, modes[n], lambda)
        denom = _trapz_inner(x, modes[n], modes[n], lambda)
        coeffs[n] = abs(denom) < 1e-12 ? 0.0 : numer / denom
    end
    return coeffs
end

"""Project normalized observations into a Gegenbauer spectral representation."""
function project_to_gegenbauer(obs::ObservationTable; n_modes::Int = 12, lambda::Float64 = 0.75)
    n_modes > 0 || throw(ArgumentError("n_modes must be positive."))
    haskey(obs.columns, :z) || throw(ArgumentError("ObservationTable must contain :z for projection."))

    z = obs.columns[:z]
    length(z) >= 3 || throw(ArgumentError("Need at least 3 vertical levels for Gegenbauer projection."))

    projected = Dict{Symbol, Vector{Float64}}()
    for key in (:u, :v, :theta, :q)
        if haskey(obs.columns, key)
            projected[key] = _project_profile_to_basis(z, obs.columns[key]; n_modes = n_modes, lambda = lambda)
        end
    end

    return (
        coefficients = projected,
        lambda = lambda,
        n_modes = n_modes,
        z_span = extrema(z),
        variables = collect(keys(projected)),
    )
end
