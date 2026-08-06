# AtmosphericSlowManifold.jl: Project Overview

`AtmosphericSlowManifold.jl` is a Julia framework for data-driven operator discovery, Geometric Singular Perturbation Theory (GSPT) analysis, spectral discretization, and parameter calibration for atmospheric boundary layer (ABL) slow manifolds and turbulence closures.

The platform establishes an end-to-end scientific pipeline: ingesting observational field data, learning symbolic operator closures subject to physical constraints, analyzing slow-manifold invariant geometry, solving prognostic atmospheric PDEs using specialized spectral Galerkin backends, and quantifying parameter uncertainty.

---

## Scientific & Execution Workflow

```
[ Tower / CSV Data ] ──► Observation ──► Manifold (State Coordinates)
                                              │
  ┌───────────────────────────────────────────┴───────────────────────────────────────────┐
  ▼                                                                                       ▼
Geometry (GSPT Analysis)                                                      Discovery (Weak WSINDy)
  • Critical Manifolds / Folds                                                  • Feature Library & Constraints
  • Canards & Desingularization                                                 • STRidge / Constrained QP
  • Continuation Kernels                                                        • Unified discover() Entrypoint
  └───────────────────────────────────────────┬───────────────────────────────────────────┘
                                              ▼
                                    Closures (WSINDyClosure)
                                              │
                                              ▼
                                     System (PDE Assembly)
                                              │
                                              ▼
                               Discretization (Spectral / FD)
                                 • Gegenbauer 3-Tensors C_ijk
                                 • Separable Modal Operators
                                              │
                                              ▼
                                  Calibration (MLE / MCMC / VI)
                                 • Point Estimates & UQ Chains

```

---

## Subsystem Architecture

### 1. Data Ingestion & Observation (`src/Observation/`)

* **Data Ingestion:** Reads NetCDF, CSV, and Cabauw tower profiles into structured `ObservationTable` objects.
* **Spectral Transformations:** Computes Gauss-Gegenbauer quadrature grid representations and spatial derivatives for boundary-layer profiles.
* **Validation:** Filters missing entries, verifies spatial monotonicity, and standardizes vertical coordinate grids $z \in [z_0, H]$.

#### Sibling Project Data Resolution (Pkgdir-Anchored)

Use package-root anchored resolution for robust cross-repository ingestion:

```julia
using AtmosphericSlowManifold

data_dir = resolve_sibling_data_dir(
  sibling_project = "SpectralBL-Analytics",
  data_subdir = "data",
)

files = find_data_files(data_dir; extensions = [".csv", ".nc"], recursive = true)
```

Load mapped columns into canonical `ObservationTable` keys (`:z`, `:u`, `:v`, `:theta`, `:q`, `:u_star`):

```julia
obs = read_observation_data(
  joinpath(data_dir, "boundary_layer_profile_01.csv");
  z_col = :height,
  u_col = :u_velocity,
  v_col = :v_velocity,
  temp_col = :potential_temp,
)

println(keys(obs.columns))
```

### 2. Manifold State Coordinates (`src/Manifold/`)

* **State Transformations:** Maps raw physical variables $(u, v, \theta, q)$ to intrinsic fast-slow manifold coordinates.
* **Dimension Reduction:** Standardizes boundary layer state vectors for geometric nullcline identification and slow-mode extraction.

### 3. Geometry Engine (`src/Geometry/`)

* **GSPT Analysis:** Analyzes fast-slow dynamics in atmospheric boundary layer equations via Geometric Singular Perturbation Theory.
* **Compiled Kernels:** Compiles symbolic `JacobianModel` evaluations via `build_function` for high-performance numerical continuation.
* **Invariant Sets:** Implements an `AbstractInvariantSet` taxonomy covering `CriticalManifoldSurface`, `FoldCurve`, `CanardSegment`, `Fenichel` normal hyperbolicity, and desingularized flow fields.

### 4. Weak Operator Discovery (`src/Discovery/`)

* **Parametric IR:** Defines symbolic `OperatorTerm{T}` building blocks and `FeatureLibrary` containers.
* **Weak-Form Assembly:** Integrates test function families (Gegenbauer test functions) over noise-corrupted observations to construct linear systems $G \boldsymbol{\Xi} \approx \mathbf{b}$ without numerical differentiation.
* **Physical Constraints:** Enforces convex inequality constraints ($A_{\text{ineq}} \boldsymbol{\Xi} \ge \mathbf{b}_{\text{ineq}}$) for realizability and monotonicity.
* **Optimization:** Solves sparse identification problems via Sequential Thresholded Ridge Regression (`STRidge`) or Quadratic Programming (`ConstrainedQP`).
* **Unified Entrypoint:** Exposes `discover(obs, library, constraints, test_family, optimizer)` to output `DiscoveredModel{Float64}` objects.

### 5. Symbolic Closures (`src/Closures/`)

* **Closure Representations:** Converts `DiscoveredModel` term sets into executable `WSINDyClosure` or standard Monin-Obukhov Similarity Theory (`MOSTClosure`) expressions.
* **Symbolic Coupling:** Binds discovered eddy diffusivities $K_m(z, u, \theta)$ directly to prognostic atmospheric PDE solvers.

### 6. Spectral & Spatial Discretization (`src/Discretization/`)

* **Finite Differences:** `MethodOfLinesFD` provides baseline spatial discretization on non-uniform stretched grids.
* **Spectral Gegenbauer Backend:** `SpectralBLGalerkin` projects boundary layer profiles onto Gegenbauer polynomials $C_n^{(\lambda)}(z)$ with orthogonal weight $w_\lambda(z) = (1 - z^2)^{\lambda - 1/2}$.
* **3-Tensor Projections:** Precomputes $C_{ijk}^{(\lambda)}$ (triple product), $A_{ijk}^{(\lambda)}$ (advection), and $B_{ijk}^{(\lambda)}$ (variable diffusivity flux) via Gauss-Gegenbauer quadrature.
* **Modal Operator Decomposition:** Separates modal RHS tendencies into explicit linear diffusion, nonlinear advection, and state-dependent diffusivity blocks with response scaling parameters (`advection_response_scale`, `diffusivity_response_scale`).
* **Diagnostics:** Exposes `ModalBudgetDiagnostic` and `evaluate_modal_budget` for per-mode transport budget decomposition.

### 7. System Assembly (`src/System/`)

* **PDE Builder:** Assembles prognostic boundary layer equations (`PrognosticPDE`) coupling momentum and scalar transport.
* **Integration Interface:** Orchestrates time-stepping via DifferentialEquations.jl with unified backend dispatches (`solve_scm`).

### 8. Calibration & UQ (`src/Calibration/`)

* **Abstract Interface:** Defines `AbstractCalibrationAlgorithm` with unified `calibrate(model, obs; algorithm=...)` dispatch.
* **Maximum Likelihood (`MaximumLikelihood`):** Closed-form MAP / Ridge estimation for fast point-estimation benchmarks and $R^2$ diagnostics.
* **Bayesian MCMC (`BayesianMCMC`):** Linear-Gaussian posterior sampling engine yielding parameter sample vectors (`:c_1`, `:c_2`, ...) for uncertainty quantification.
* **Variational Inference (`VariationalInference`):** Mean-field Gaussian VI backend with ELBO tracing, convergence diagnostics, and posterior draw generation.

### 9. Multi-Format Export Utilities (`src/System/ExportUtilities.jl`)

* **CSV Export:** Writes `DataFrame`, `ObservationTable`, and `ModalBudgetDiagnostic` outputs.
* **JSON Export:** Serializes discovered operator models, sparse coefficients, and calibration diagnostics.
* **NetCDF Export:** Writes 2D $(z,t)$ trajectory fields with dimension metadata for downstream analysis.

---

## Subsystem Implementation Status

### Markdown Summary

| Subsystem Module | Primary Capabilities | Implementation Status | Test Coverage |
| --- | --- | --- | --- |
| **`Observation`** | NetCDF/CSV ingestion, Gegenbauer transforms | **Complete** | Pass |
| **`Manifold`** | Intrinsic state representation, coordinate transforms | **Complete** | Pass |
| **`Geometry`** | Fast-slow GSPT, compiled Jacobian kernels, fold/canard geometry | **Complete** | Pass |
| **`Discovery`** | Weak WSINDy, feature libraries, constraints, unified `discover()` | **Complete** | Pass |
| **`Closures`** | `WSINDyClosure`, `MOSTClosure`, symbolic PDE coupling | **Complete** | Pass |
| **`Discretization`** | Stretched FD & Gegenbauer Galerkin backends, 3-tensors $C_{ijk}^{(\lambda)}$ | **Complete** | Pass |
| **`System`** | Boundary layer PDE assembly, time-stepping orchestration | **Complete** | Pass |
| **`Calibration`** | Unified interface, `MaximumLikelihood`, `BayesianMCMC`, `VI` | **Complete** | Pass |
| **`System/ExportUtilities`** | CSV/JSON/NetCDF export for model and trajectory artifacts | **Complete** | Pass |

### LaTeX Summary

```latex
\begin{table}[htbp]
  \centering
  \caption{Subsystem Implementation Status and Test Verification}
  \label{tab:subsystem-status}
  \begin{tabular}{llcc}
    \toprule
    \textbf{Subsystem Module} & \textbf{Primary Capabilities} & \textbf{Status} & \textbf{Test Suite} \\
    \midrule
    \textbf{\texttt{Observation}} & NetCDF/CSV ingestion, Gegenbauer transforms & Complete & Pass \\
    \textbf{\texttt{Manifold}} & Intrinsic state representation, coordinates & Complete & Pass \\
    \textbf{\texttt{Geometry}} & Fast--slow GSPT, Jacobian kernels, folds/canards & Complete & Pass \\
    \textbf{\texttt{Discovery}} & Weak WSINDy, constraints, \texttt{discover()} & Complete & Pass \\
    \textbf{\texttt{Closures}} & \texttt{WSINDyClosure}, \texttt{MOSTClosure}, PDE coupling & Complete & Pass \\
    \textbf{\texttt{Discretization}} & Stretched FD \& Gegenbauer Galerkin, 3-tensors & Complete & Pass \\
    \textbf{\texttt{System}} & Boundary layer PDE assembly, time-stepping & Complete & Pass \\
    	extbf{\texttt{Calibration}} & Unified interface, \texttt{MaximumLikelihood}, \texttt{BayesianMCMC}, \texttt{VariationalInference} & Complete & Pass \\
    	extbf{\texttt{System/ExportUtilities}} & CSV/JSON/NetCDF artifact export utilities & Complete & Pass \\
    \bottomrule
  \end{tabular}
\end{table}

```

---

## Verification & Test Architecture

The package test suite covers unit, module, and integration-level behavior across all 8 subsystems:

1. **`test_wsindy_engine.jl` & `test_discovery_modules.jl`:** Asserts weak-form assembly, sparse regression solvers (`STRidge`, `ConstrainedQP`), and unified `discover()` execution paths.
2. **`test_scm_backends.jl`:** Verifies finite-difference and spectral Galerkin solvers, 3-tensor symmetries ($C_{kij} = C_{kji}$), modal operator block decompositions, and solver-level trajectory divergence under `ASM_RUN_SMOKE=1`.
3. **`test_geometry_*.jl`:** Validates GSPT invariant manifolds, fold curve desingularization, and compiled Jacobian performance.
4. **`test_calibration_interface.jl`:** Asserts deterministic coefficient recovery ($R^2 \approx 1.0$), posterior sample vectors (`BayesianMCMC`), and ELBO-backed variational outputs (`VariationalInference`).
5. **`test_export_utilities.jl`:** Verifies artifact export correctness across CSV, JSON, and NetCDF pathways.