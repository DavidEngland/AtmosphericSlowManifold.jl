module EnergyBudget

using Statistics

export SurfaceEnergySummary, surface_energy_budget, tke_budget, closure_dissipation, energy_residual

struct SurfaceEnergySummary{T<:Real}
    mean_imbalance::T
    rms_imbalance::T
    max_imbalance::T
    cumulative_imbalance::T
end

"""
    surface_energy_budget(R_n, H, G, LE)

Evaluates surface energy closure residual over time:
`imbalance = R_n - H - G - LE`.
"""
function surface_energy_budget(
    R_n::AbstractVector{T},
    H::AbstractVector{T},
    G::AbstractVector{T},
    LE::AbstractVector{T},
) where {T<:Real}
    n = length(R_n)
    (n == length(H) && n == length(G) && n == length(LE)) || throw(ArgumentError("Array dimensions must match."))
    n == 0 && return SurfaceEnergySummary{Float64}(0.0, 0.0, 0.0, 0.0)

    acc = 0.0
    acc2 = 0.0
    cum = 0.0
    max_abs = 0.0
    @inbounds for i in eachindex(R_n, H, G, LE)
        imb = Float64(R_n[i]) - Float64(H[i]) - Float64(G[i]) - Float64(LE[i])
        aimb = abs(imb)
        acc += imb
        acc2 += imb * imb
        cum += imb
        if aimb > max_abs
            max_abs = aimb
        end
    end

    return SurfaceEnergySummary{Float64}(
        acc / n,
        sqrt(acc2 / n),
        max_abs,
        cum,
    )
end

"""
    tke_budget(u, v)

Computes pointwise turbulent kinetic energy proxy `0.5*(u^2 + v^2)`.
"""
function tke_budget(u::AbstractVector{T}, v::AbstractVector{T}) where {T<:Real}
    length(u) == length(v) || throw(ArgumentError("Array dimensions must match."))
    out = similar(u, Float64)
    @inbounds for i in eachindex(u, v)
        ui = Float64(u[i])
        vi = Float64(v[i])
        out[i] = 0.5 * (ui * ui + vi * vi)
    end
    return out
end

"""
    closure_dissipation(K_m, du_dz)

Calculates local turbulent production/dissipation surrogate:
`epsilon_m(z) = K_m(z) * (du_dz)^2`.
"""
function closure_dissipation(K_m::AbstractVector{T}, du_dz::AbstractVector{T}) where {T<:Real}
    length(K_m) == length(du_dz) || throw(ArgumentError("Array dimensions must match."))
    out = similar(K_m, Float64)
    @inbounds for i in eachindex(K_m, du_dz)
        g = Float64(du_dz[i])
        out[i] = Float64(K_m[i]) * g * g
    end
    return out
end

"""
    energy_residual(u, v, z_grid)

Computes trapezoid-integrated kinetic energy column
`integral 0.5*(u^2 + v^2) dz`.
"""
function energy_residual(u::AbstractVector{T}, v::AbstractVector{T}, z_grid::AbstractVector{T}) where {T<:Real}
    n = length(z_grid)
    (n == length(u) && n == length(v)) || throw(ArgumentError("Array dimensions must match."))
    n <= 1 && return 0.0

    integrated_ke = 0.0
    @inbounds for i in 1:(n - 1)
        u1 = Float64(u[i])
        v1 = Float64(v[i])
        u2 = Float64(u[i + 1])
        v2 = Float64(v[i + 1])
        ke1 = 0.5 * (u1 * u1 + v1 * v1)
        ke2 = 0.5 * (u2 * u2 + v2 * v2)
        dz = Float64(z_grid[i + 1]) - Float64(z_grid[i])
        integrated_ke += 0.5 * (ke1 + ke2) * dz
    end
    return integrated_ke
end

end # module EnergyBudget
