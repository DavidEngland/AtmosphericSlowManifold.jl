# src/Closures/MOSTClosure.jl
"""
    MOSTClosure{T<:AbstractFloat} <: AbstractClosure

Monin-Obukhov Similarity Theory baseline closure parameterization.

Default coefficients follow modern SBL practice:
- kappa = 0.40
- pr_t = 1.0
- cm_stable = 4.7

Classical Businger-Dyer style unstable forms are included through
`gamma_m` and `gamma_h` defaults (16.0, 16.0).
"""
struct MOSTClosure{T<:AbstractFloat} <: AbstractClosure
    kappa::T
    pr_t::T
    cm_stable::T
    gamma_m::T
    gamma_h::T
end

function MOSTClosure(
    ;
    kappa::Real = 0.40,
    pr_t::Real = 1.0,
    cm_stable::Real = 4.7,
    gamma_m::Real = 16.0,
    gamma_h::Real = 16.0,
)
    T = promote_type(typeof(kappa), typeof(pr_t), typeof(cm_stable), typeof(gamma_m), typeof(gamma_h))
    Tf = T <: AbstractFloat ? T : Float64
    return MOSTClosure{Tf}(Tf(kappa), Tf(pr_t), Tf(cm_stable), Tf(gamma_m), Tf(gamma_h))
end

MOSTClosure(kappa::Real, pr_t::Real, cm_stable::Real) =
    MOSTClosure(; kappa = Float64(kappa), pr_t = Float64(pr_t), cm_stable = Float64(cm_stable))

MOSTClosure(kappa::Real, pr_t::Real, cm_stable::Real, gamma_m::Real, gamma_h::Real) =
    MOSTClosure(
        ;
        kappa = Float64(kappa),
        pr_t = Float64(pr_t),
        cm_stable = Float64(cm_stable),
        gamma_m = Float64(gamma_m),
        gamma_h = Float64(gamma_h),
    )

@inline function _most_local_height(m::ManifoldState)
    return hasproperty(m, :z) ? m.z : m.z0
end

@inline function _most_phi_m(c::MOSTClosure{T}, r) where {T<:AbstractFloat}
    stable = one(T) + c.cm_stable * r
    unstable = (one(T) - c.gamma_m * r)^(-T(0.25))
    return ifelse(r >= zero(T), stable, unstable)
end

@inline function _most_phi_h(c::MOSTClosure{T}, r) where {T<:AbstractFloat}
    stable = c.pr_t * (one(T) + c.cm_stable * r)
    unstable = c.pr_t * (one(T) - c.gamma_h * r)^(-T(0.5))
    return ifelse(r >= zero(T), stable, unstable)
end

function eddy_momentum(c::MOSTClosure{T}, m::ManifoldState) where {T<:AbstractFloat}
    z_eval = _most_local_height(m)
    phi_m = max(_most_phi_m(c, m.r), T(0.1))
    return (c.kappa * m.u_star * z_eval) / phi_m
end

function eddy_heat(c::MOSTClosure{T}, m::ManifoldState) where {T<:AbstractFloat}
    z_eval = _most_local_height(m)
    phi_h = max(_most_phi_h(c, m.r), T(0.1))
    return (c.kappa * m.u_star * z_eval) / phi_h
end

surface_flux(::MOSTClosure, m::ManifoldState) = m.u_star^2

"""
    evaluate_diffusivity_profile!(K_out, c, z_grid; r_profile=nothing, u_star=0.3, z0=0.0, mask_below_z0=false)

In-place evaluation of MOST momentum diffusivity over `z_grid`.

Example
```julia
K = zeros(Float64, 64)
z = collect(range(0.0, 100.0; length = 64))
c = MOSTClosure()
evaluate_diffusivity_profile!(K, c, z; r_profile = fill(0.05, 64), u_star = 0.3)
```
"""
function evaluate_diffusivity_profile!(
    K_out::AbstractVector{T},
    c::MOSTClosure{T},
    z_grid::AbstractVector{T},
    ;
    r_profile::Union{Nothing, AbstractVector{T}} = nothing,
    u_star::T = T(0.3),
    z0::T = zero(T),
    mask_below_z0::Bool = false,
) where {T<:AbstractFloat}
    length(K_out) == length(z_grid) || throw(DimensionMismatch("K_out and z_grid must have equal length."))
    r_profile === nothing || (length(r_profile) == length(z_grid) || throw(DimensionMismatch("r_profile and z_grid must have equal length.")))

    @inbounds @simd for i in eachindex(K_out, z_grid)
        z_raw = z_grid[i]
        if mask_below_z0 && z_raw <= z0
            K_out[i] = zero(T)
        else
            z = max(z_raw, z0)
            r = r_profile === nothing ? zero(T) : r_profile[i]
            phi_m = max(_most_phi_m(c, r), T(0.1))
            K_out[i] = (c.kappa * u_star * z) / phi_m
        end
    end
    return K_out
end

function evaluate_heat_diffusivity_profile!(
    K_out::AbstractVector{T},
    c::MOSTClosure{T},
    z_grid::AbstractVector{T},
    ;
    r_profile::Union{Nothing, AbstractVector{T}} = nothing,
    u_star::T = T(0.3),
    z0::T = zero(T),
    mask_below_z0::Bool = false,
) where {T<:AbstractFloat}
    length(K_out) == length(z_grid) || throw(DimensionMismatch("K_out and z_grid must have equal length."))
    r_profile === nothing || (length(r_profile) == length(z_grid) || throw(DimensionMismatch("r_profile and z_grid must have equal length.")))

    @inbounds @simd for i in eachindex(K_out, z_grid)
        z_raw = z_grid[i]
        if mask_below_z0 && z_raw <= z0
            K_out[i] = zero(T)
        else
            z = max(z_raw, z0)
            r = r_profile === nothing ? zero(T) : r_profile[i]
            phi_h = max(_most_phi_h(c, r), T(0.1))
            K_out[i] = (c.kappa * u_star * z) / phi_h
        end
    end
    return K_out
end
