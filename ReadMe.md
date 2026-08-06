# AtmosphericSlowManifold.jl

`AtmosphericSlowManifold.jl` is a reusable Julia scientific computing framework designed for turbulence closure discovery, Geometric Singular Perturbation Theory (GSPT) analysis, and reduced-order prognostic modeling of fast–slow geophysical continuum systems—specifically the Atmospheric Boundary Layer (ABL).

The framework decouples observation ingestion, manifold state representation, GSPT geometric algorithms, sparse operator identification, symbolic PDE assembly, and spatial discretization into independent, orthogonal layers.

---

## Architectural Principles

```
  [ Observation Space ] (Tower Data, LES, Radiosondes, NetCDF)
            │
            ▼
  [ Intrinsic Manifold ] (η₁, η₂, η₃, R, Ω, χ, Π_G, λ_min)
            │
            ▼
  [ Geometric Engine ] ──► (Critical Manifolds, Folds, Canards, Fenichel, Continuation)
            │
            ▼
  [ Symbolic Discovery ] ─► (Weak-form WSINDy, Parametric IR OperatorTerm{T}, Constraints)
            │
            ▼
    [ Prognostic PDE ] ───► (ModelingToolkit.jl System Assembly)
            │
            ▼
 [ Spatial Discretization ] (Method-of-Lines FD, Spectral Gegenbauer Galerkin)

```

1. **Orthogonality of Spaces:** Physical state coordinates, geometric manifold analysis, symbolic closure mathematical structure, prognostic PDEs, and spatial discretizations are fully decoupled.
2. **Autonomous Geometric Engine:** The `Geometry` layer computes invariant manifolds, fold curves, and normal hyperbolicity persistence bounds without depending on the discovery engine or numerical solvers.
3. **Discretization-Independent Operator Discovery:** Turbulence closures are identified as symbolic operators in continuous space via weak-form integration, independent of downstream numerical grids.
4. **Parametric Intermediate Representation:** Discovered operators are represented as typed `OperatorTerm{T}` structures ($T \in \{\text{Float64}, \text{Num}, \text{Dual}, \text{Measurement}\}$) before compilation into `ModelingToolkit.jl` Abstract Syntax Trees (ASTs).

---

## Subsystem Overview

```
AtmosphericSlowManifold.jl/
├── src/
│   ├── Observation/           # Data ingestion & spectral transforms (Cabauw, NetCDF)
│   ├── Manifold/              # Intrinsic manifold coordinates & state representations
│   ├── Geometry/              # Autonomous GSPT algorithms & continuation
│   │   ├── Jacobians.jl       # Symbolics build_function compiled numerical kernels
│   │   ├── CriticalManifold.jl# Fast nullcline surface root-solver (0 = f(x, y, p))
│   │   ├── DesingularizedFlow.jl # Rescaled slow dynamics & singular point scanners
│   │   ├── Fenichel.jl        # Spectral gap & normal hyperbolicity metrics
│   │   ├── FoldTracking.jl    # Bordered Jacobian continuation for det(J_y) = 0
│   │   ├── CanardDetection.jl # Folded singularity classification (nodes, saddles, foci)
│   │   └── Continuation.jl    # Abstract continuation API & pseudo-arclength loops
│   ├── Discovery/             # Weak-form sparse operator discovery (WSINDy)
│   │   ├── SymbolicExtraction.jl # OperatorTerm{T} IR & MTK AST generator
│   │   ├── LibraryBuilder.jl  # Typed feature library (State, Derivative, Diagnostic)
│   │   ├── ConstraintBuilder.jl # Physical inequality builders (Positivity, Monotonicity)
│   │   ├── TestFunctions.jl   # Gegenbauer, B-Spline space-time test bases
│   │   ├── WeakForms.jl       # Quadrature & integration-by-parts (G Ξ ≈ b)
│   │   ├── SparseRegression.jl# STRidge & JuMP-backed convex QP solvers
│   │   └── WSINDyEngine.jl    # Discovery pipeline orchestrator
│   ├── Closures/              # Symbolic closures (WSINDyClosure, MOSTClosure)
│   │   ├── PrognosticPDE.jl   # ModelingToolkit boundary-layer PDE assembly
│   │   └── SurfaceBoundary.jl # Surface energy balance & flux boundary conditions
│   ├── Discretization/        # Interchangeable numerical backends
│   │   ├── Backends/MethodOfLinesFD.jl # Stretched finite differences
│   │   └── Backends/SpectralBLGalerkin.jl # Boundary layer Gegenbauer Galerkin
│   └── Calibration/           # Hierarchical Turing.jl MCMC inference
└── test/                      # Comprehensive test suite across all subsystems

```

---

## Installation

`AtmosphericSlowManifold.jl` requires Julia 1.9 or higher.

```julia
using Pkg
Pkg.add(url="https://github.com/davidengland/AtmosphericSlowManifold.jl.git")

```

Or clone locally and instantiate:

```bash
git clone https://github.com/davidengland/AtmosphericSlowManifold.jl.git
cd AtmosphericSlowManifold.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'

```

---

## Quickstart

### 1. Geometric Singular Perturbation Theory (GSPT) Analysis

Analyze fast–slow dynamics, compile numerical Jacobian models, and track critical manifolds:

```julia
using AtmosphericSlowManifold
using AtmosphericSlowManifold.Geometry
using ModelingToolkit, Symbolics

# Define fast-slow boundary layer system
@variables u(t) v(t) w(t) Ri(t)
@parameters Ug Pi_G

# System of fast algebraic nullclines and slow equations
f_fast = [Ug - u + w*v, -v - w*u]
fast_vars = [u, v]
slow_vars = [w, Ri]
params = [Ug, Pi_G]

# Compile JacobianModel with Symbolics build_function kernels
jmodel = JacobianModel(f_fast, fast_vars, slow_vars, params)

# Compute critical manifold surface where f_fast = 0
cm_surface = compute_critical_manifold(jmodel, [10.0, 0.0], [0.1, 0.2], [10.0, 0.01])

# Calculate Fenichel normal hyperbolicity metrics across the manifold
fenichel_diag = profile_fenichel_hyperbolicity(cm_surface, jmodel, [10.0, 0.01])
println("Min Spectral Gap: ", fenichel_diag.min_spectral_gap)

```

### 2. Discovered Turbulence Closure Integration

Build space-time weak forms from observation data, assemble physical constraints, run sparse regression, and execute the PDE system:

```julia
using AtmosphericSlowManifold
using AtmosphericSlowManifold.Observation
using AtmosphericSlowManifold.Discovery
using AtmosphericSlowManifold.System
using AtmosphericSlowManifold.Discretization

# 1. Ingest observation data
obs = read_tower_csv("data/cabauw_profile.csv")

# 2. Build feature library and weak-form matrices
lib = build_feature_library([:u, :theta], [:Pi_G], 2)
test_basis = GegenbauerFamily(0.5, 10)
weak_sys = assemble_weak_system(obs, test_basis, lib)

# 3. Enforce physical constraints (e.g., eddy diffusivity K_m >= 0)
positivity = PositivityConstraint(StateVariable(:u))
constraints = assemble_constraint_matrix([positivity], lib, obs.evaluation_grid)

# 4. Sparse identification using JuMP-backed Quadratic Programming
coeffs = solve_sparse_regression(weak_sys.G, weak_sys.b,
                                 constraints.A_ineq, constraints.b_ineq,
                                 ConstrainedQP(lambda=1e-4))

# 5. Extract parametric OperatorTerm IR and generate ModelingToolkit AST
discovered_model = DiscoveredModel(:K_m, coeffs, lib)
closure = build_wsindy_closure(discovered_model)

# 6. Assemble prognostic PDE and solve via Stretched Finite Differences
pde_system = build_pde_system(closure)
backend = MethodOfLinesFD(n_points=50, stretch_factor=2.5)
solution = solve_scm(pde_system, closure, backend, (0.0, 86400.0))

```

---

## Testing

Run the test suite to verify all subsystem contracts:

```julia
using Pkg
Pkg.test("AtmosphericSlowManifold")

```

The test suite covers:

* `test_geometry_core.jl` & `test_geometry_foundations.jl`: Critical manifolds, fold curves, canards, continuation algorithms, and compiled Jacobian models.
* `test_discovery_ir.jl`: Typed `OperatorTerm{T}` IR conversions, feature libraries, and physical constraint assemblies.
* `test_wsindy_engine.jl`: Weak-form matrix assembly and sparse recovery algorithms.
* `test_scm_backends.jl`: Spatial discretization backends and boundary layer prognostic PDE solutions.
* `test_observation_ingestion.jl` & `test_spectral_transform.jl`: Data pipeline and Gegenbauer transforms.

---

## Production Campaign Exports

The repository includes a production batch export pipeline for campaign artifacts derived from sibling `SpectralBL-Analytics` data sources.

### Command Reference

```bash
# Execute end-to-end extraction across CASES-99, FLOSS, BLLAST, and SHEBA
make campaign-export

# Display the generated multi-campaign summary table in the terminal
make campaign-summary

# Validate generated campaign artifacts and report schemas
make campaign-validate

# Remove generated campaign artifact directories
make campaign-clean

```

Generate all current outputs for CASES-99, FLOSS, BLLAST, and SHEBA:

```bash
make campaign-export

```

This runs [scripts/run_campaign_exports.jl](scripts/run_campaign_exports.jl) and writes artifacts under [reports/generated/campaign_exports](reports/generated/campaign_exports).

Available artifact groups:

* [reports/generated/campaign_exports/csv](reports/generated/campaign_exports/csv): raw campaign extracts and compact diagnostic summary tables
* [reports/generated/campaign_exports/json](reports/generated/campaign_exports/json): discovered-model-style coefficient payloads and serialized diagnostics
* [reports/generated/campaign_exports/netcdf](reports/generated/campaign_exports/netcdf): exported 2D stability / Richardson-number fields
* [reports/generated/campaign_exports/figures](reports/generated/campaign_exports/figures): campaign metric plots, Ri heatmaps where available, and a comparative campaign overview figure
* [reports/generated/campaign_exports/tables/campaign_overview.csv](reports/generated/campaign_exports/tables/campaign_overview.csv): derived per-campaign overview metrics including observation count, inferred height levels, mean wind speed, and mean stability
* [reports/generated/campaign_exports/tables/campaign_summary.md](reports/generated/campaign_exports/tables/campaign_summary.md): markdown manifest of generated outputs
* [reports/generated/campaign_exports/tables/campaign_summary.tex](reports/generated/campaign_exports/tables/campaign_summary.tex): LaTeX report table

To print the current generated manifest directly:

```bash
make campaign-summary

```

### Target Overview

| Target | Executable Command | Primary Output / Action |
|---|---|---|
| `make campaign-export` | `julia --project=. scripts/run_campaign_exports.jl` | Exports CSV, JSON, NetCDF, and PNG artifacts under `reports/generated/campaign_exports/` |
| `make campaign-summary` | `cat reports/generated/campaign_exports/tables/campaign_summary.md` | Prints campaign stats and generated manifest summary |
| `make campaign-validate` | `julia --project=. scripts/validate_campaign_exports.jl` | Verifies required artifacts, overview schema, and report sections |
| `make campaign-clean` | `rm -rf reports/generated/campaign_exports/` | Removes generated campaign outputs |

---

## License

This project is licensed under the MIT License.