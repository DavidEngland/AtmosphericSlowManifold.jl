# src/Closures/Z0HRClosure.jl

"""
    Z0HRClosure{T<:AbstractFloat} <: AbstractClosure

Zero-Offset Hyperbolic Regularization (Z0HR) closure parameterization.
Eliminates branch-switching discontinuities across neutral (Ri = 0) and
critical (Ri = Ri_c) thresholds using C^∞ hyperbolic coordinates.

Default parameters:
- kappa = 0.40
- pr_t = 1.0
- B_um = 16.0
- B_uh = 16.0
- Ri_c = 0.25
- eps = 1e-3
"""
struct Z0HRClosure{T<:AbstractFloat} <: AbstractClosure
    kappa::T
    pr_t::T
    B_um::T
    B_uh::T
    Ri_c::T
    eps::T
end

function Z0HRClosure(;
    kappa::Real=0.40,
    pr_t::Real=1.0,
    B_um::Real=16.0,
    B_uh::Real=16.0,
    Ri_c::Real=0.25,
    eps::Real=1e-3,
)
    T = promote_type(typeof(kappa), typeof(pr_t), typeof(B_um), typeof(B_uh), typeof(Ri_c), typeof(eps))
    Tf = T <: AbstractFloat ? T : Float64
    return Z0HRClosure{Tf}(Tf(kappa), Tf(pr_t), Tf(B_um), Tf(B_uh), Tf(Ri_c), Tf(eps))
end

@inline function _z0hr_local_height(m::ManifoldState)
    return hasproperty(m, :z) ? m.z : m.z0
end

@inline function z0hr_stability_functions(c::Z0HRClosure{T}, Ri) where {T<:AbstractFloat}
    Ri_pos = smooth_max(Ri; eps = c.eps)
    Ri_neg = smooth_min(Ri; eps = c.eps)

    g_raw = one(T) - (Ri_pos / c.Ri_c)
    stable_factor = smooth_max(g_raw; eps = c.eps)

    unstab_m = one(T) - (c.B_um * Ri_neg)
    unstab_h = one(T) - (c.B_uh * Ri_neg)

    S_m = (stable_factor^2) * sqrt(unstab_m)
    S_h = (one(T) / c.pr_t) * (stable_factor^2) * (unstab_h^T(0.75))

    return S_m, S_h
end

function eddy_momentum(c::Z0HRClosure{T}, m::ManifoldState) where {T<:AbstractFloat}
    z_eval = _z0hr_local_height(m)
    S_m, _ = z0hr_stability_functions(c, m.r)
    return c.kappa * m.u_star * z_eval * S_m
end

function eddy_heat(c::Z0HRClosure{T}, m::ManifoldState) where {T<:AbstractFloat}
    z_eval = _z0hr_local_height(m)
    _, S_h = z0hr_stability_functions(c, m.r)
    return c.kappa * m.u_star * z_eval * S_h
end

surface_flux(::Z0HRClosure, m::ManifoldState) = m.u_star^2

function evaluate_diffusivity_profile!(
    K_out::AbstractVector{T},
    c::Z0HRClosure{T},
    z_grid::AbstractVector{T};
    r_profile::Union{Nothing,AbstractVector{T}}=nothing,
    u_star::T=T(0.3),
    z0::T=zero(T),
    mask_below_z0::Bool=false,
) where {T<:AbstractFloat}
    length(K_out) == length(z_grid) || throw(DimensionMismatch("K_out and z_grid must have equal length."))
    r_profile === nothing || (length(r_profile) == length(z_grid) || throw(DimensionMismatch("r_profile and z_grid must have equal length.")))

    @inbounds @simd for i in eachindex(K_out, z_grid)
        z_raw = z_grid[i]
        if mask_below_z0 && z_raw <= z0
            K_out[i] = zero(T)
        else
            z = smooth_floor(z_raw, z0; eps=c.eps)
            r = r_profile === nothing ? zero(T) : r_profile[i]
            S_m, _ = z0hr_stability_functions(c, r)
            K_out[i] = c.kappa * u_star * z * S_m
        end
    end
    return K_out
end

function evaluate_heat_diffusivity_profile!(
    K_out::AbstractVector{T},
    c::Z0HRClosure{T},
    z_grid::AbstractVector{T};
    r_profile::Union{Nothing,AbstractVector{T}}=nothing,
    u_star::T=T(0.3),
    z0::T=zero(T),
    mask_below_z0::Bool=false,
) where {T<:AbstractFloat}
    length(K_out) == length(z_grid) || throw(DimensionMismatch("K_out and z_grid must have equal length."))
    r_profile === nothing || (length(r_profile) == length(z_grid) || throw(DimensionMismatch("r_profile and z_grid must have equal length.")))

    @inbounds @simd for i in eachindex(K_out, z_grid)
        z_raw = z_grid[i]
        if mask_below_z0 && z_raw <= z0
            K_out[i] = zero(T)
        else
            z = smooth_floor(z_raw, z0; eps=c.eps)
            r = r_profile === nothing ? zero(T) : r_profile[i]
            _, S_h = z0hr_stability_functions(c, r)
            K_out[i] = c.kappa * u_star * z * S_h
        end
    end
    return K_out
end