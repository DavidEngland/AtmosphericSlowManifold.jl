The provided framework document synthesizes the mathematical foundation and software architecture of `AtmosphericSlowManifold.jl` (ASM.jl). All seven core functional domains mapped in the specification correspond directly to the implemented modules in `src/`, exports in `AtmosphericSlowManifold.jl`, and workflow targets in the `Makefile`.

---

### System Architecture & Implementation Mapping Matrix

| Architectural Layer | Source Module (`src/`) | Key Types & Functions | Computational & Physical Role |
| --- | --- | --- | --- |
| **1. Manifold & Geometry** | `src/Geometry/` | `ManifoldState`, `FoldConstraint`, `compute_jacobian` | Tracks critical manifold $\mathcal{S}_0$, calculates fast-subsystem Jacobian $D_{\mathbf{y}}\mathbf{g}$, detects fold loci ($\det D_{\mathbf{y}}\mathbf{g} = 0$), and tracks Fenichel normal hyperbolicity loss. |
| **2. Observation & Data Ingestion** | `src/Observation/` | `ObservationTable`, `read_tower_csv`, `project_to_gegenbauer` | Ingests campaign NetCDF/CSV datasets (CASES-99, SHEBA, FLOSS, BLLAST), resolves variable aliases ($H, u_*$), and projects raw profiles onto Gegenbauer bases. |
| **3. Symbolic Discovery Engine** | `src/Discovery/` | `WSINDyEngine`, `ConstraintBuilder`, `ModelSelection` | Assembles weak-form space-time projections $(\mathbf{G}, \mathbf{b})$, enforces physical constraints ($K_m, K_h \ge 0$, energy monotonicity), solves STRidge/JuMP QPs, and evaluates AIC/BIC/Pareto fronts. |
| **4. Turbulence Closures** | `src/Closures/` | `AbstractAtmosphericClosure`, `MOSTClosure`, `WSINDyClosure` | Houses $C^\infty$-differentiable closure formulations. Uses smooth floors and hyperbolic weights to eliminate derivative step-discontinuities at $Ri_{\text{cr}}$. |
| **5. Discretization & Solvers** | `src/Discretization/` | `MethodOfLinesFD`, `SpectralBLGalerkin`, `generate_stretched_grid` | Offers stretched 1D finite differences and zero-allocation Gegenbauer-Galerkin spectral expansions using precomputed 3-tensor contraction matrices ($C_{ijk}^{(\lambda)}$). |
| **6. System & Calibration** | `src/Calibration/` | `MaximumLikelihood`, `BayesianMCMC`, `VariationalInference`, `HierarchicalTuring`, `UncertaintyPropagation` | Converts discovered symbols to ModelingToolkit.jl systems. Conducts MLE point estimation, Turing.jl MCMC sampling, mean-field VI, multi-site hierarchical calibration, and profile credibility ribbon evaluations. |
| **7. Diagnostics & Reporting** | `src/Diagnostics/`<br>

<br>`src/Discovery/LaTeXExporter.jl` | `ErrorMetrics`, `EnergyBudget`, `to_latex`, `latex_term_table` | Aggregates statistical metrics ($R^2$, RMSE), verifies energy budgets, and auto-compiles discovered PDEs and summary statistics into LaTeX tables and equations. |

---

### Core Mathematical Mechanics Integrated in `src/`

#### 1. GSPT Adjugate Rescaling & Canard Trajectory Resolution

Near fold lines where $\det(D_{\mathbf{y}}\mathbf{f}) = 0$, standard flow vectors diverge. The `src/Geometry/` engine rescales the independent time variable ($d\tau = dt / \det(D_{\mathbf{y}}\mathbf{f})$) and applies Cramer's rule using the adjugate matrix:

$$\frac{d\mathbf{z}}{d\tau} = \operatorname{adj}(D_{\mathbf{y}}\mathbf{f}) \cdot \mathbf{g}(\mathbf{z})$$

This desingularization enables smooth integration through folded singularities (saddles, nodes, foci) and classifies non-equilibrium canard trajectories crossing from stable to unstable sheets.

#### 2. Weak-Form WSINDy Space-Time Projection

Noise robustness is guaranteed by projecting the governing PDE $\mathcal{L}[u] = 0$ against smooth, compactly supported Gegenbauer or B-spline test functions $\phi_{i,j}(z, t)$:

$$\int_{\Omega_t} \int_{\Omega_z} \mathcal{L}[u](z, t) \, \phi_{i,j}(z, t) \, dz \, dt = 0$$

Integrating by parts transfers spatial and temporal derivatives from noisy observational data $u(z, t)$ onto $C^\infty$ test functions $\phi_{i,j}$, constructing the linear system $\mathbf{G}\boldsymbol{\xi} = \mathbf{b}$ without numerical differentiation.

#### 3. $C^\infty$ Smooth Operators for Continuity & Persistence

To prevent singular Jacobians and preserve Fenichel persistence during continuation solver steps, step discontinuities are replaced with smooth algebraic floor functions:

$$\operatorname{smooth\_max}(x, y; \epsilon) = \frac{x + y + \sqrt{(x - y)^2 + \epsilon^2}}{2}$$

---

### Pipeline Execution Architecture

```
                       make pde-benchmark / make campaign-export
                                           │
         ┌─────────────────────────────────┼─────────────────────────────────┐
         ▼                                 ▼                                 ▼
[ 1. Discovery & Selection ]     [ 2. Bayesian Calibration ]      [ 3. Automated Reporting ]
 • Weak-form projections (G, b)   • MLE / Turing MCMC / VI         • LaTeX equations (to_latex)
 • Positivity & Energy QP bounds  • Hierarchical multi-site pooling  • Term & summary tables (.tex)
 • AIC / BIC & Pareto extraction  • 95% K_m(z) profile ribbons     • Publication JSON summaries

```

The underlying library in `src/` is complete and tested. Enhancing `scripts/run_pde_closure_benchmark.jl` is the next step to output these mathematical models into compiled LaTeX tables and Pareto metrics under `reports/`.