module Diagnostics

include("ErrorMetrics.jl")
include("EnergyBudget.jl")
include("ManifoldMetrics.jl")
include("CampaignDiagnostics.jl")

using .ErrorMetrics
using .EnergyBudget
using .ManifoldMetrics
using .CampaignDiagnosticsModule

export ErrorMetrics, EnergyBudget, ManifoldMetrics, CampaignDiagnosticsModule
export CampaignDiagnostics
export rmse, mae, bias, r2, nrmse, skill_score, correlation,
       normalized_bias, relative_l2_error, closure_residual
export SurfaceEnergySummary, surface_energy_budget, tke_budget, closure_dissipation, energy_residual
export transversality, fold_distance, slow_manifold_error, normal_hyperbolicity

end # module Diagnostics
