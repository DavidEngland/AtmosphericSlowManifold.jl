<!-- Auto-generated from package source -->
> **Source:** `src/Diagnostics/CampaignDiagnostics.md`

`CampaignDiagnosticsModule` defines the `CampaignDiagnostics` parametric struct, providing a unified container that aggregates statistical error metrics, thermodynamic energy balances, and geometric slow-manifold diagnostics for evaluation against observational field campaign datasets (e.g., CASES-99, GABLS1).

---

### Campaign Diagnostics Aggregation Schema

| Field Name | Type | Source Diagnostic Submodule | Metric Description & Physical Context |
| --- | --- | --- | --- |
| **`campaign_name`** | `String` | Metadata | Campaign identifier or station dataset tag. |
| **`num_samples`** | `Int` | Metadata | Total vertical profile or temporal sample count. |
| **`rmse_u`** | `T` | `ErrorMetrics` | Root Mean Square Error for horizontal velocity profile $u(z)$. |
| **`rmse_theta`** | `T` | `ErrorMetrics` | Root Mean Square Error for potential temperature profile $\theta(z)$. |
| **`r2_closure`** | `T` | `ErrorMetrics` | Coefficient of determination $R^2$ for learned or baseline closure fits. |
| **`mean_closure_residual`** | `T` | `ErrorMetrics` | Pointwise mean absolute residual between predicted and observed eddy fluxes. |
| **`energy_imbalance_rms`** | `T` | `EnergyBudget` | RMS surface energy balance residual ($R_n - H - G - LE$) across time series. |
| **`integrated_dissipation`** | `T` | `EnergyBudget` | Vertically integrated column turbulent dissipation rate ($\int \varepsilon_m \, dz$). |
| **`transversality_mean`** | `T` | `ManifoldMetrics` | Mean alignment angle $\tau$ between fast tendencies and slow-manifold normals. |
| **`transversality_std`** | `T` | `ManifoldMetrics` | Standard deviation of fast-slow transversality across samples. |
| **`manifold_distance_mean`** | `T` | `ManifoldMetrics` | Mean Euclidean distance $\text{SME}$ between observed profiles and slow manifold projections. |
| **`fold_distance_min`** | `T` | `ManifoldMetrics` | Minimum distance to nearest fold locus, indicating proximity to stability transition limits. |

---

### End-to-End Aggregation Workflow

```julia
using AtmosphericSlowManifold.Diagnostics

# 1. Compute individual metrics from profile inputs
val_rmse_u   = ErrorMetrics.rmse(u_pred, u_obs)
val_rmse_th  = ErrorMetrics.rmse(theta_pred, theta_obs)
e_summary    = EnergyBudget.surface_energy_budget(Rn, H, G, LE)
tau_val      = ManifoldMetrics.transversality(v_fast, n_slow)
dist_fold    = ManifoldMetrics.fold_distance(state, fold_locus)

# 2. Package into unified CampaignDiagnostics summary
summary = CampaignDiagnostics(
    "CASES-99_Night3",
    length(u_obs),
    val_rmse_u,
    val_rmse_th,
    0.942,
    0.015,
    e_summary.rms_imbalance,
    4.12,
    tau_val,
    0.03,
    0.18,
    dist_fold
)

```

This structural consolidation ensures that evaluation reports exported via `ExportUtilities.export_to_json` capture both raw statistical performance and underlying manifold stability diagnostics in a single record.