## Diagnostics Subsystem Architecture (src/Diagnostics/)  
The frozen numerical evaluation framework decouples state integration from quantitative validation, physical energy balance, and slow-manifold geometry.  
```
src/Diagnostics/
├── Diagnostics.jl             # Submodule re-exports & unifying API
├── ErrorMetrics.jl            # Predictive accuracy & skill scores
├── EnergyBudget.jl            # Conservation, surface flux balance, & dissipation
├── ManifoldMetrics.jl         # Geometric transversality, hyperbolicity, & fold distance
└── CampaignDiagnostics.jl     # Unified multi-campaign comparison structure

```
**1. src/Diagnostics/ErrorMetrics.jl**  
```
module ErrorMetrics

using Statistics
using LinearAlgebra

export rmse, mae, bias, r2, nrmse, skill_score, correlation, 
       normalized_bias, relative_l2_error, closure_residual

"""
    rmse(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where T<:Real

Computes Root Mean Square Error:
\$\\text{RMSE} = \\sqrt{\\frac{1}{N} \\sum_{i=1}^N (y_{\\text{pred},i} - y_{\\text{true},i})^2}\$
"""
function rmse(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where T<:Real
    @assert length(y_pred) == length(y_true) "Array dimensions must match"
    return sqrt(mean(abs2, y_pred .- y_true))
end

"""
    mae(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where T<:Real

Computes Mean Absolute Error:
\$\\text{MAE} = \\frac{1}{N} \\sum_{i=1}^N |y_{\\text{pred},i} - y_{\\text{true},i}|\$
"""
function mae(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where T<:Real
    @assert length(y_pred) == length(y_true) "Array dimensions must match"
    return mean(abs, y_pred .- y_true)
end

"""
    bias(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where T<:Real

Computes Mean Systematic Bias:
\$\\text{Bias} = \\frac{1}{N} \\sum_{i=1}^N (y_{\\text{pred},i} - y_{\\text{true},i})\$
"""
function bias(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where T<:Real
    @assert length(y_pred) == length(y_true) "Array dimensions must match"
    return mean(y_pred .- y_true)
end

"""
    r2(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where T<:Real

Computes Coefficient of Determination (\$R^2\$):
\$R^2 = 1 - \\frac{\\sum (y_{\\text{true},i} - y_{\\text{pred},i})^2}{\\sum (y_{\\text{true},i} - \\bar{y}_{\\text{true}})^2}\$
"""
function r2(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where T<:Real
    @assert length(y_pred) == length(y_true) "Array dimensions must match"
    ss_res = sum(abs2, y_true .- y_pred)
    y_bar = mean(y_true)
    ss_tot = sum(abs2, y_true .- y_bar)
    return iszero(ss_tot) ? one(T) : one(T) - (ss_res / ss_tot)
end

"""
    nrmse(y_pred::AbstractArray{T}, y_true::AbstractArray{T}; norm::Symbol=:std) where T<:Real

Computes Normalized RMSE scaled by `:std`, `:range`, or `:mean`.
"""
function nrmse(y_pred::AbstractArray{T}, y_true::AbstractArray{T}; norm::Symbol=:std) where T<:Real
    val_rmse = rmse(y_pred, y_true)
    scale = if norm == :std
        std(y_true)
    elseif norm == :range
        maximum(y_true) - minimum(y_true)
    elseif norm == :mean
        abs(mean(y_true))
    else
        error("Unsupported normalization option: $norm")
    end
    return iszero(scale) ? val_rmse : val_rmse / scale
end

"""
    skill_score(y_pred::AbstractArray{T}, y_true::AbstractArray{T}, y_ref::AbstractArray{T}) where T<:Real

Computes Murphy skill score relative to a baseline reference trajectory (e.g. unclosed neutral model):
\$\\text{Skill} = 1 - \\frac{\\text{MSE}(\\text{pred}, \\text{true})}{\\text{MSE}(\\text{ref}, \\text{true})}\$
"""
function skill_score(y_pred::AbstractArray{T}, y_true::AbstractArray{T}, y_ref::AbstractArray{T}) where T<:Real
    mse_pred = mean(abs2, y_pred .- y_true)
    mse_ref  = mean(abs2, y_ref .- y_true)
    return iszero(mse_ref) ? zero(T) : one(T) - (mse_pred / mse_ref)
end

"""
    relative_l2_error(y_pred::AbstractVector{T}, y_true::AbstractVector{T}) where T<:Real

Computes Relative \$L_2\$ Error Norm:
\$e_{L_2} = \\frac{\\|y_{\\text{pred}} - y_{\\text{true}}\\|_2}{\\|y_{\\text{true}}\\|_2}\$
"""
function relative_l2_error(y_pred::AbstractVector{T}, y_true::AbstractVector{T}) where T<:Real
    norm_true = norm(y_true)
    return iszero(norm_true) ? norm(y_pred - y_true) : norm(y_pred - y_true) / norm_true
end

"""
    closure_residual(phi_obs::AbstractVector{T}, phi_model::AbstractVector{T}) where T<:Real

Computes point-by-point closure residual vector:
\$R_k = |\\phi_{\\text{obs}}(\\zeta_k) - \\phi_{\\text{model}}(\\zeta_k)|\$
"""
function closure_residual(phi_obs::AbstractVector{T}, phi_model::AbstractVector{T}) where T<:Real
    @assert length(phi_obs) == length(phi_model) "Dimensions must match"
    return abs.(phi_obs .- phi_model)
end

end # module ErrorMetrics

```
**2. src/Diagnostics/EnergyBudget.jl**  
```
module EnergyBudget

using Statistics

export surface_energy_budget, tke_budget, closure_dissipation, energy_residual

struct SurfaceEnergySummary{T<:Real}
    mean_imbalance::T
    rms_imbalance::T
    max_imbalance::T
    cumulative_imbalance::T
end

"""
    surface_energy_budget(R_n::AbstractVector{T}, H::AbstractVector{T}, 
                         G::AbstractVector{T}, LE::AbstractVector{T}) where T<:Real

Evaluates surface energy closure residual over time:
\$\\text{Imbalance} = R_n - H - G - LE\$
"""
function surface_energy_budget(R_n::AbstractVector{T}, H::AbstractVector{T}, 
                               G::AbstractVector{T}, LE::AbstractVector{T}) where T<:Real
    imbalance = R_n .- H .- G .- LE
    return SurfaceEnergySummary(
        mean(imbalance),
        sqrt(mean(abs2, imbalance)),
        maximum(abs, imbalance),
        sum(imbalance)
    )
end

"""
    closure_dissipation(K_m::AbstractVector{T}, du_dz::AbstractVector{T}) where T<:Real

Calculates local turbulent energy production/dissipation rate:
\$\\varepsilon_m(z) = K_m(z) \\left( \\frac{\\partial u}{\\partial z} \\right)^2\$
"""
function closure_dissipation(K_m::AbstractVector{T}, du_dz::AbstractVector{T}) where T<:Real
    return K_m .* (du_dz .^ 2)
end

"""
    energy_residual(u::AbstractVector{T}, v::AbstractVector{T}, 
                    z_grid::AbstractVector{T}, dt::T) where T<:Real

Computes column-integrated domain kinetic energy time derivative \$\\frac{dE}{dt}\$.
"""
function energy_residual(u::AbstractVector{T}, v::AbstractVector{T}, 
                         z_grid::AbstractVector{T}) where T<:Real
    # Trapezoidal integration of 0.5 * (u^2 + v^2)
    ke_density = @. 0.5 * (u^2 + v^2)
    integrated_ke = zero(T)
    @inbounds for i in 1:(length(z_grid)-1)
        dz = z_grid[i+1] - z_grid[i]
        integrated_ke += 0.5 * (ke_density[i] + ke_density[i+1]) * dz
    end
    return integrated_ke
end

end # module EnergyBudget

```
**3. src/Diagnostics/ManifoldMetrics.jl**  
```
module ManifoldMetrics

using LinearAlgebra
using Statistics

export transversality, fold_distance, slow_manifold_error, normal_hyperbolicity

"""
    transversality(v_fast::AbstractVector{T}, n_slow::AbstractVector{T}) where T<:Real

Calculates alignment between fast velocity tendencies and slow-manifold unit normals:
\$\\mathcal{T} = \\frac{|v_{\\text{fast}} \\cdot n_{\\text{slow}}|}{\\|v_{\\text{fast}}\\|_2 \\|n_{\\text{slow}}\\|_2}\$
"""
function transversality(v_fast::AbstractVector{T}, n_slow::AbstractVector{T}) where T<:Real
    norm_product = norm(v_fast) * norm(n_slow)
    return iszero(norm_product) ? zero(T) : abs(dot(v_fast, n_slow)) / norm_product
end

"""
    fold_distance(state::AbstractVector{T}, fold_locus::AbstractMatrix{T}) where T<:Real

Computes Euclidean distance from a trajectory point to the nearest fold singularity locus.
"""
function fold_distance(state::AbstractVector{T}, fold_locus::AbstractMatrix{T}) where T<:Real
    min_dist = typemax(T)
    @inbounds for i in 1:size(fold_locus, 2)
        dist = norm(state .- @view(fold_locus[:, i]))
        if dist < min_dist
            min_dist = dist
        end
    end
    return min_dist
end

"""
    slow_manifold_error(state::AbstractVector{T}, manifold_map::Function) where T<:Real

Evaluates residual distance between prognostic state and mapped slow invariant manifold manifold_map(state).
"""
function slow_manifold_error(state::AbstractVector{T}, manifold_map::Function) where T<:Real
    return norm(state .- manifold_map(state))
end

"""
    normal_hyperbolicity(jacobian::AbstractMatrix{T}, n_fast::Int) where T<:Real

Computes the spectral gap ratio between fast and slow Jacobian eigenvalues:
\$\\Gamma = \\frac{\\min |\\text{Re}(\\lambda_{\\text{fast}})|}{\\max |\\text{Re}(\\lambda_{\\text{slow}})|}\$
"""
function normal_hyperbolicity(jacobian::AbstractMatrix{T}, n_fast::Int) where T<:Real
    evals = sort(abs.(real.(eigen(jacobian).values)), rev=true)
    fast_rate = evals[n_fast]
    slow_rate = evals[n_fast + 1]
    return iszero(slow_rate) ? typemax(T) : fast_rate / slow_rate
end

end # module ManifoldMetrics

```
**4. src/Diagnostics/CampaignDiagnostics.jl**  
```
module CampaignDiagnosticsModule

using ..ErrorMetrics
using ..EnergyBudget
using ..ManifoldMetrics

export CampaignDiagnostics

struct CampaignDiagnostics{T<:AbstractFloat}
    campaign_name::String
    num_samples::Int
    # Statistical Metrics
    rmse_u::T
    rmse_theta::T
    r2_closure::T
    mean_closure_residual::T
    # Physical / Conservation
    energy_imbalance_rms::T
    integrated_dissipation::T
    # Geometric Slow-Manifold Metrics
    transversality_mean::T
    transversality_std::T
    manifold_distance_mean::T
    fold_distance_min::T
end

end # module CampaignDiagnosticsModule

```
**Single VS Code Implementation Prompt**  
```
Task: Implement the Diagnostics subsystem in AtmosphericSlowManifold.jl following the exact 4-module layout.

Locations:
- Create `src/Diagnostics/ErrorMetrics.jl`
- Create `src/Diagnostics/EnergyBudget.jl`
- Create `src/Diagnostics/ManifoldMetrics.jl`
- Create `src/Diagnostics/CampaignDiagnostics.jl`
- Create `src/Diagnostics/Diagnostics.jl` (includes all 4 modules and re-exports top-level symbols)
- Include `src/Diagnostics/Diagnostics.jl` inside `src/AtmosphericSlowManifold.jl` and re-export `Diagnostics`.
- Create unit tests in `test/test_diagnostics.jl` and include in `test/runtests.jl`.

Requirements:
1. Ensure full type stability (no `Any`, zero runtime allocations inside inner metric loops).
2. Use `@inbounds` and `@views` for all array slice operations.
3. Write test cases covering RMSE/MAE/R2, surface energy balance calculations, and transversality functions.
4. Verify that `julia --project=. -e 'using AtmosphericSlowManifold; using Test; include("test/test_diagnostics.jl")'` passes with 100% green assertions.

```
