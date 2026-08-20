"""
    PhysicalSimilarityClosure{T<:AbstractFloat} <: AbstractClosure

Physical closure model parameterized by discovered similarity diagnostics.

Fields
- `phi_coeffs`: polynomial coefficients for `phi_m(zeta)` in ascending order.
- `zeta_coeffs`: polynomial coefficients for transformed stability parameter in ascending order.
- `karman`: Von Karman constant.
- `ustar`: friction velocity magnitude.
- `L_obukhov`: Obukhov length scale.
- `z_ref`: reference height used in campaign parameterization.
"""
struct PhysicalSimilarityClosure{T<:AbstractFloat} <: AbstractClosure
    phi_coeffs::Vector{T}
    zeta_coeffs::Vector{T}
    karman::T
    ustar::T
    L_obukhov::T
    z_ref::T
end

@inline function _poly_eval(coeffs::AbstractVector{T}, x::Number) where {T<:AbstractFloat}
    isempty(coeffs) && return one(promote_type(T, typeof(x)))
    acc = coeffs[end]
    @inbounds for i in (length(coeffs)-1):-1:1
        acc = muladd(acc, x, coeffs[i])
    end
    return acc
end

@inline function _zeta_from_height(c::PhysicalSimilarityClosure{T}, z::Number) where {T<:AbstractFloat}
    raw = z / c.L_obukhov
    transformed = isempty(c.zeta_coeffs) ? raw : _poly_eval(c.zeta_coeffs, raw)
    # Smoothly bound zeta to keep phi(zeta) evaluations well-conditioned in extreme regimes.
    zeta_limit = T(5.0)
    return zeta_limit * tanh(transformed / zeta_limit)
end

@inline function _phi_m(c::PhysicalSimilarityClosure{T}, zeta::Number) where {T<:AbstractFloat}
    return smooth_floor(_poly_eval(c.phi_coeffs, zeta), T(0.1); eps = T(1e-3))
end

# Minimum physical background diffusivity so strongly stable regimes (Ri_g >> 0.25)
# cannot drive K_m/K_h to zero and stall implicit time-stepping.
const K_MIN_DIFFUSIVITY = 1.0e-4

@inline function _km(c::PhysicalSimilarityClosure{T}, z::Number) where {T<:AbstractFloat}
    zeta = _zeta_from_height(c, z)
    return max(T(K_MIN_DIFFUSIVITY), (c.karman * c.ustar * z) / _phi_m(c, zeta))
end

"""
    PhysicalSimilarityClosure(; phi_coeffs, zeta_coeffs, karman=0.4, ustar=0.3, L_obukhov=-100.0, z_ref=10.0)

Construct a physical similarity closure directly from numeric parameters.
"""
function PhysicalSimilarityClosure(
    ;
    phi_coeffs::AbstractVector{<:Real} = [1.0],
    zeta_coeffs::AbstractVector{<:Real} = [0.0, 1.0],
    karman::Real = 0.4,
    ustar::Real = 0.3,
    L_obukhov::Real = -100.0,
    z_ref::Real = 10.0,
)
    T = Float64
    L = isfinite(Float64(L_obukhov)) && abs(Float64(L_obukhov)) > eps(Float64) ? Float64(L_obukhov) : -100.0
    return PhysicalSimilarityClosure{T}(
        Float64.(collect(phi_coeffs)),
        Float64.(collect(zeta_coeffs)),
        Float64(karman),
        Float64(ustar),
        L,
        Float64(z_ref),
    )
end

@inline function _json_number(x, default::Float64)
    if x === nothing
        return default
    elseif x isa Number
        v = Float64(x)
        return isfinite(v) ? v : default
    else
        s = strip(String(x))
        isempty(s) && return default
        v = tryparse(Float64, s)
        return (v === nothing || !isfinite(v)) ? default : v
    end
end

function _json_get(obj, keys::Tuple{Vararg{String}}, default)
    cur = obj
    for k in keys
        if cur === nothing
            return default
        end
        try
            if haskey(cur, k)
                cur = cur[k]
            else
                return default
            end
        catch
            return default
        end
    end
    return cur
end

function _term_coefficient(payload, term_name::String, default::Float64)
    terms = _json_get(payload, ("terms",), nothing)
    terms === nothing && return default
    for term in terms
        name = _json_get(term, ("name",), nothing)
        name === nothing && continue
        String(name) == term_name || continue
        return _json_number(_json_get(term, ("coefficient",), nothing), default)
    end
    return default
end

function _resolve_diagnostics_json_path(json_path::String)
    if isfile(json_path)
        return json_path
    end

    base = basename(json_path)
    diagnostics_dir = joinpath(pwd(), "reports", "generated", "campaign_exports", "json")

    candidates = String[]
    push!(candidates, joinpath(diagnostics_dir, base))

    if occursin("_diagnostics.json", base)
        alt = replace(base, "_diagnostics.json" => "_model_and_diagnostics.json")
        push!(candidates, joinpath(dirname(json_path), alt))
        push!(candidates, joinpath(diagnostics_dir, alt))
    elseif occursin("_model_and_diagnostics.json", base)
        alt = replace(base, "_model_and_diagnostics.json" => "_diagnostics.json")
        push!(candidates, joinpath(dirname(json_path), alt))
        push!(candidates, joinpath(diagnostics_dir, alt))
    end

    for c in candidates
        isfile(c) && return c
    end
    throw(ArgumentError("Could not find diagnostics JSON at $(json_path) or compatible fallback paths."))
end

"""
    PhysicalSimilarityClosure(json_path::String)

Construct a physical similarity closure from a campaign JSON diagnostics payload.
The parser supports current generated artifacts under
`reports/generated/campaign_exports/json/*_model_and_diagnostics.json`.
"""
function PhysicalSimilarityClosure(json_path::String)
    resolved_path = _resolve_diagnostics_json_path(json_path)
    payload = JSON3.read(read(resolved_path, String))

    ustar = _json_number(_json_get(payload, ("diagnostics", "stats", "ustar", "mean"), nothing), 0.3)
    L_mean = _json_number(_json_get(payload, ("diagnostics", "obukhov_scaling", "mean"), nothing), -100.0)
    L_mean = abs(L_mean) > eps(Float64) ? L_mean : -100.0

    phi_mean = _json_number(_json_get(payload, ("diagnostics", "similarity_parameters", "phi_obs", "mean"), nothing), 1.0)
    zeta_mean = _json_number(_json_get(payload, ("diagnostics", "similarity_parameters", "zeta", "mean"), nothing), 0.1)

    phi_term = _term_coefficient(payload, "phi_obs", 0.0)
    zeta_term = _term_coefficient(payload, "zeta", 1.0)

    phi0 = max(phi_mean, 0.1)
    phi1 = isfinite(phi_term) ? phi_term : 0.0
    zref_from_stats = _json_number(_json_get(payload, ("diagnostics", "stats", "z", "mean"), nothing), NaN)
    zref = if isfinite(zref_from_stats)
        zref_from_stats
    else
        zeta_abs = max(abs(zeta_mean), 1e-6)
        abs(L_mean) * zeta_abs
    end

    return PhysicalSimilarityClosure(
        phi_coeffs = [phi0, phi1],
        zeta_coeffs = [zeta_mean, zeta_term],
        karman = 0.4,
        ustar = ustar,
        L_obukhov = L_mean,
        z_ref = zref,
    )
end

"""
    (closure::PhysicalSimilarityClosure)(z::Real)

Evaluate momentum eddy diffusivity at height `z`.
"""
@inline function (closure::PhysicalSimilarityClosure{T})(z::Real) where {T<:AbstractFloat}
    zT = T(z)
    z_eff = smooth_floor(zT, zero(T); eps = T(1e-3))
    return T(_km(closure, z_eff))
end

"""
    evaluate_diffusivity_profile!(K_out, closure, z_grid)

In-place evaluation of momentum eddy diffusivity over a vertical grid.
"""
function evaluate_diffusivity_profile!(
    K_out::AbstractVector{T},
    closure::PhysicalSimilarityClosure{T},
    z_grid::AbstractVector{T},
) where {T<:AbstractFloat}
    length(K_out) == length(z_grid) || throw(DimensionMismatch("K_out and z_grid must have equal length."))
    @views @inbounds for i in eachindex(z_grid)
        K_out[i] = closure(z_grid[i])
    end
    return K_out
end

"""
    evaluate_heat_diffusivity_profile!(K_out, closure, z_grid)

In-place evaluation of heat eddy diffusivity over a vertical grid.
"""
function evaluate_heat_diffusivity_profile!(
    K_out::AbstractVector{T},
    closure::PhysicalSimilarityClosure{T},
    z_grid::AbstractVector{T},
) where {T<:AbstractFloat}
    length(K_out) == length(z_grid) || throw(DimensionMismatch("K_out and z_grid must have equal length."))
    @views @inbounds for i in eachindex(z_grid)
        z = smooth_floor(z_grid[i], zero(T); eps = T(1e-3))
        zeta = _zeta_from_height(closure, z)
        phi_h = smooth_floor(_phi_m(closure, zeta) * T(1.1), T(0.1); eps = T(1e-3))
        K_out[i] = max(T(K_MIN_DIFFUSIVITY), (closure.karman * closure.ustar * z) / phi_h)
    end
    return K_out
end

function eddy_momentum(c::PhysicalSimilarityClosure, m::ManifoldState)
    z_eval = hasproperty(m, :z) ? m.z : m.z0
    return _km(c, z_eval)
end

function eddy_heat(c::PhysicalSimilarityClosure, m::ManifoldState)
    z_eval = hasproperty(m, :z) ? m.z : m.z0
    return _km(c, z_eval) / 1.1
end

surface_flux(c::PhysicalSimilarityClosure, ::ManifoldState) = c.ustar^2
