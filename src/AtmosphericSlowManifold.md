The `AtmosphericSlowManifold.jl` root module serves as the orchestration hub for the framework, explicitly mapping the modular subcomponents into a unified API surface while relying on core Julia ecosystem tools (`ModelingToolkit`, `JuMP`, `MethodOfLines`, `DifferentialEquations`).

```
AtmosphericSlowManifold.jl (Root Module)
 ├── Manifold/ & Geometry/   --> ManifoldState, FoldConstraint, Geometry
 ├── Observation/            --> ObservationTable, SpectralBLTransform (Gegenbauer)
 ├── Discovery/              --> WSINDy, FeatureLibrary, PhysicalConstraints, STRidge/ConstrainedQP
 ├── Closures/               --> AbstractClosure (WSINDy, MOST, PhysicalSimilarity)
 ├── Discretization/         --> MethodOfLinesFD, SpectralBLGalerkin, StretchedGrid
 ├── System/ & Diagnostics/  --> SurfaceBoundary, PrognosticPDE (MTK), Diagnostics
 └── Calibration/ & Export/  --> BayesianMCMC, MaximumLikelihood, VariationalInference, ExportUtilities

```

### Architectural Mapping & Module Breakdown

1. **Manifold State & Geometry (`Manifold/`, `Geometry/`)**
* **Exports:** `ManifoldState`, `FoldConstraint`, `fold_residual`, `fold_transversality`, `Geometry`.
* **Function:** Encapsulates the GSPT fast-slow state transformations and exposes diagnostic functions to compute fold proximity ($\det D_{\mathbf{y}}\mathbf{g} = 0$) and manifold transversality ($\bar{\mathcal{T}}$).


2. **Observational Data Ingestion (`Observation/`)**
* **Exports:** `ObservationTable`, `read_tower_csv`, `read_tower_netcdf`, `project_to_gegenbauer`.
* **Function:** Normalizes tower and field campaign data (CASES-99, FLOSS, SHEBA, BLLAST) into standardized table layouts and projects vertical profiles onto modal Gegenbauer bases.


3. **Weak-Form Symbolic Discovery Engine (`Discovery/`)**
* **Exports:** `build_weak_library`, `fit_wsindy_jump`, `GegenbauerFamily`, `BSplineFamily`, `STRidge`, `ConstrainedQP`, `PositivityConstraint`, `EnergyConstraint`, `assemble_constraint_matrix`.
* **Function:** Implements weak-form WSINDy regression. Integrates `JuMP.jl` for hard-constrained quadratic programming (`ConstrainedQP`) to enforce positivity ($K_m, K_h \ge 0$), energy dissipation, and asymptotic bounds during feature selection.


4. **Closure Abstractions (`Closures/`)**
* **Exports:** `AbstractAtmosphericClosure`, `WSINDyClosure`, `MOSTClosure`, `PhysicalSimilarityClosure`, `eddy_momentum`, `eddy_heat`, `evaluate_diffusivity_profile!`.
* **Function:** Provides a standard API for diffusivity profiles and surface fluxes, allowing empirical baselines (MOST) and learned symbolic closures (WSINDy) to be interchanged.


5. **Discretization & Numerical Backends (`Discretization/`)**
* **Exports:** `MethodOfLinesFD`, `SpectralBLGalerkin`, `generate_stretched_grid`, `SpectralNonlinearTensors`, `precompute_nonlinear_tensors`, `solve_scm`.
* **Function:** Decouples the spatial discretization from the physics. Supports grid-stretching finite differences via `MethodOfLines.jl` or zero-allocation modal projections via Gegenbauer contraction tensors (`SpectralNonlinearTensors`).


6. **Prognostic PDE Assembly & Calibration (`System/`, `Calibration/`, `ExportUtilities/`)**
* **Exports:** `build_pde_system`, `to_mtk_expression`, `BayesianMCMC`, `MaximumLikelihood`, `VariationalInference`, `calibrate`.
* **Function:** Translates discovered symbolic expressions (`to_mtk_expression`) directly into `ModelingToolkit.jl` prognostic PDE systems (`build_pde_system`). Handles parameter estimation and serialization to CSV, JSON, or NetCDF.



---

### Key System Integration Highlights

* **Symbolic Interoperability:** Discovered operators in `Discovery` export directly to `ModelingToolkit.jl` expressions (`to_mtk_expression`), enabling automated generation of Jacobian matrices and stiff ODE integration via standard `DifferentialEquations.jl` solvers.
* **Optimized Execution Paths:** Pre-computing non-linear tensors (`precompute_nonlinear_tensors`) in `Discretization/Backends/SpectralBLGalerkin.jl` allows modal RHS evaluations to run without dynamic memory allocations in time-stepping loops.
* **Unified Interface Layer:** Solvers (`solve_scm`) and calibration dispatchers (`dispatch_calibrate`) operate uniformly regardless of whether a finite-difference or spectral Galerkin backend is specified.

Would you like to examine the inner implementation of a specific component, such as `Discovery/Discovery.jl` for the constrained sparse regression assembly, or `Discretization/Backends/SpectralBLGalerkin.jl` for the tensor contraction algebra?