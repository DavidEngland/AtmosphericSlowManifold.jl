# src/Closures/MOSTClosure.jl
"""
    MOSTClosure{T<:AbstractFloat} <: AbstractClosure

Monin-Obukhov Similarity Theory baseline closure parameterization.
"""
struct MOSTClosure{T<:AbstractFloat} <: AbstractClosure
    kappa::T
    pr_t::T
    cm_stable::T
end

function MOSTClosure(; kappa::Real = 0.41, pr_t::Real = 0.74, cm_stable::Real = 4.7)
    T = promote_type(typeof(kappa), typeof(pr_t), typeof(cm_stable))
    Tf = T <: AbstractFloat ? T : Float64
    return MOSTClosure{Tf}(Tf(kappa), Tf(pr_t), Tf(cm_stable))
end

MOSTClosure(kappa::Real, pr_t::Real, cm_stable::Real) =
    MOSTClosure(; kappa = Float64(kappa), pr_t = Float64(pr_t), cm_stable = Float64(cm_stable))

@inline function _most_local_height(m::ManifoldState)
    return hasproperty(m, :z) ? m.z : m.z0
end

function eddy_momentum(c::MOSTClosure{T}, m::ManifoldState) where {T<:AbstractFloat}
    z_eval = _most_local_height(m)
    stability_factor = max(one(T) + c.cm_stable * m.r, T(0.1))
    return (c.kappa * m.u_star * z_eval) / stability_factor
end

eddy_heat(c::MOSTClosure{T}, m::ManifoldState) where {T<:AbstractFloat} = eddy_momentum(c, m) / c.pr_t
surface_flux(::MOSTClosure, m::ManifoldState) = m.u_star^2

function evaluate_diffusivity_profile!(
    K_out::AbstractVector{T},
    c::MOSTClosure{T},
    z_grid::AbstractVector{T},
) where {T<:AbstractFloat}
    length(K_out) == length(z_grid) || throw(DimensionMismatch("K_out and z_grid must have equal length."))
    # Neutral-to-stable MOST baseline in profile form with r=0 surrogate.
    denom = max(one(T), T(0.1))
    @inbounds for i in eachindex(K_out, z_grid)
        z = max(z_grid[i], zero(T))
        K_out[i] = (c.kappa * T(0.3) * z) / denom
    end
    return K_out
end

function evaluate_heat_diffusivity_profile!(
    K_out::AbstractVector{T},
    c::MOSTClosure{T},
    z_grid::AbstractVector{T},
) where {T<:AbstractFloat}
    evaluate_diffusivity_profile!(K_out, c, z_grid)
    @inbounds for i in eachindex(K_out)
        K_out[i] /= c.pr_t
    end
    return K_out
end
