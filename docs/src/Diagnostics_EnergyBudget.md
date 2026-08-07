<!-- Auto-generated from package source -->
> **Source:** `src/Diagnostics/EnergyBudget.md`

The `EnergyBudget` module provides numerical diagnostics for evaluating thermodynamic and mechanical energy conservation across boundary layer simulations, including surface energy balance, kinetic energy profiles, and eddy dissipation rates.

---

### Functional Overview & Formulations

| Diagnostic Function | Mathematical Formulation | Return Type | Physical Interpretation |
| --- | --- | --- | --- |
| **`surface_energy_budget`** | $\text{imb}_i = R_{n,i} - H_i - G_i - LE_i$ | `SurfaceEnergySummary` | Evaluates surface energy closure imbalance (net radiation $R_n$, sensible heat flux $H$, ground flux $G$, latent heat flux $LE$) and returns mean, RMS, maximum, and cumulative imbalances. |
| **`tke_budget`** | $e_i = \frac{1}{2}(u_i^2 + v_i^2)$ | `Vector{Float64}` | Calculates local horizontal kinetic energy density proxy vector along grid profile points. |
| **`closure_dissipation`** | $\varepsilon_{m,i} = K_{m,i} \left(\frac{\partial u}{\partial z}\right)_i^2$ | `Vector{Float64}` | Computes local turbulent momentum dissipation surrogate from eddy diffusivity $K_m$ and vertical shear $\frac{\partial u}{\partial z}$. |
| **`energy_residual`** | $E_{\text{col}} = \int_{0}^{z_{\text{top}}} \frac{1}{2}(u^2 + v^2) \, dz$ | `Float64` | Integrates column kinetic energy density over vertical spatial grid $z$ using 1D trapezoidal integration: <br>$$\sum_{i=1}^{N-1} \frac{e_i + e_{i+1}}{2} (z_{i+1} - z_i)$$

 |

---

### Key Numerical & Implementation Features

* **Trapezoidal Column Integration:** `energy_residual` integrates non-uniform vertical grids $z$ without requiring fixed spacing $\Delta z$, allowing diagnostic tracking on stretched boundary layer meshes.
* **Array Bounds Safety:** All exported routines execute strict input length checks (`ArgumentError`) prior to entering computational loops, preventing out-of-bounds indexing or silent length mismatched evaluations.
* **Type-Safe `Float64` Accumulation:** Explicit casting inside `@inbounds` loops guarantees double-precision numerical precision when calculating cumulative energy flux sums over multi-day time series.