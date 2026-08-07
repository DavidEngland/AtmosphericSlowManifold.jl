The `Diagnostics` module acts as the central diagnostic hub in `AtmosphericSlowManifold.jl`, re-exporting evaluation metrics across four distinct submodules covering statistical verification, physical energy budgets, geometrical manifold analysis, and campaign dataset diagnostics.

---

### Diagnostic Submodule Architecture

| Submodule | Exported Types & Functions | Core Purpose |
| --- | --- | --- |
| **`ErrorMetrics`** | `rmse`, `mae`, `bias`, `r2`, `nrmse`, `skill_score`, `correlation`, `normalized_bias`, `relative_l2_error`, `closure_residual` | Statistical model validation, goodness-of-fit assessment, and residual tracking against observations or reference solutions. |
| **`EnergyBudget`** | `SurfaceEnergySummary`, `surface_energy_budget`, `tke_budget`, `closure_dissipation`, `energy_residual` | Conservation checking, surface energy balance evaluation, Turbulent Kinetic Energy (TKE) budget decomposition, and closure dissipation tracking. |
| **`ManifoldMetrics`** | `transversality`, `fold_distance`, `slow_manifold_error`, `normal_hyperbolicity` | Geometric and dynamical systems analysis tracking slow manifold distance, stability loss, fold proximity, and normal hyperbolicity breakdown. |
| **`CampaignDiagnosticsModule`** | `CampaignDiagnostics` | Aggregated multi-station observational dataset diagnostics and field campaign evaluation metrics. |

---

### Key Operational Characteristics

* **Re-exported Namespaces:** Users can access diagnostic routines directly from `Diagnostics` (e.g., `Diagnostics.rmse(...)`) or via the parent module `AtmosphericSlowManifold` without requiring explicit sub-imports.
* **Modular Separation:** Keeps statistical, thermodynamic, and differential-geometric metrics strictly decoupled while maintaining a single import point for analysis pipelines and reporting utilities.