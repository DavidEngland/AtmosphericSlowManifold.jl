`src/Geometry/Geometry.jl` forms the computational core for Geometric Singular Perturbation Theory (GSPT) in `AtmosphericSlowManifold.jl`. It encapsulates submodules that handle fast-slow subsystem decomposition, critical manifold construction, fold curve continuation, Fenichel normal hyperbolicity diagnostics, and canard orbit detection.

```
Geometry Module
├── Jacobians.jl           --> Fast Jacobian (D_y g), adjugates, tangent space bases
├── CriticalManifold.jl    --> Manifold surface root-finding (g(x,y) = 0)
├── DesingularizedFlow.jl  --> Time-rescaled slow flow vector fields
├── FoldTracking.jl        --> Fold curve tracing & indicator functions
├── CanardDetection.jl     --> Folded singularity classification & canard segments
├── Fenichel.jl            --> Normal hyperbolicity diagnostics & spectral metrics
└── Continuation.jl        --> Pseudo-arclength continuation engine

```

---

### Key Component Submodules & Mathematical Roles

| Submodule | Core Types & Methods | Mathematical / Physical Function |
| --- | --- | --- |
| **`Jacobians.jl`** | `JacobianModel`, `JacobianCache`, `compute_fast_jacobian`, `evaluate_adjugate` | Computes fast subsystem Jacobian $D_{\mathbf{y}}\mathbf{g}(\mathbf{z})$, its determinant, and adjugate matrices using direct or finite difference (`finite_difference_jacobian_y`) methods. |
| **`CriticalManifold.jl`** | `ManifoldPoint`, `CriticalManifoldSurface`, `solve_critical_surface` | Solves for the slow invariant set $\mathcal{S}_0 = \{(\mathbf{x}, \mathbf{y}) \mid \mathbf{g}(\mathbf{x}, \mathbf{y}) = 0\}$ defining equilibrium velocity and temperature profiles. |
| **`DesingularizedFlow.jl`** | `desingularized_vector_field`, `slow_flow_vector`, `find_desingularized_singular_points` | Rescales independent variable near fold points ($d\tau = dt/\det(D_{\mathbf{y}}\mathbf{g})$) to project smooth vector fields onto singular loci. |
| **`FoldTracking.jl`** | `FoldCurve`, `fold_indicator`, `track_fold_curve` | Detects zero-crossings of $\det(D_{\mathbf{y}}\mathbf{g})$ and tracks 1D fold boundaries in parameter and state space. |
| **`CanardDetection.jl`** | `FoldedSingularity`, `CanardSegment`, `classify_folded_singularity` | Identifies folded saddles, nodes, and foci, detecting canard trajectories that cross fold lines from stable to unstable manifold branches. |
| **`Fenichel.jl`** | `HyperbolicityReport`, `fenichel_metrics`, `hyperbolicity_profile` | Evaluates spectral gaps and Fenichel normal hyperbolicity metrics to bound manifold persistence under small turbulent perturbations ($\varepsilon > 0$). |
| **`Continuation.jl`** | `PseudoArclength`, `ContinuationBranch`, `continue_manifold` | Implements predictor-corrector continuation to trace multi-dimensional invariant manifolds across non-linear bifurcation boundaries. |

---

### Practical Application

This module provides the geometric engine required to diagnose boundary layer transitions. By combining `Jacobians` and `FoldTracking`, `ASM.jl` detects the onset of turbulence collapse (loss of normal hyperbolicity) without relying on arbitrary bulk Richardson number criteria ($Ri_{\text{cr}}$).