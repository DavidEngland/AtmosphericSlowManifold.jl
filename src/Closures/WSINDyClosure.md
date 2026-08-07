`src/Closures/WSINDyClosure.jl` implements the concrete data-driven closure type for `AtmosphericSlowManifold.jl`. It encapsulates symbolic mathematical expressions discovered from observational data via Weak Sparse Identification of Non-Linear Dynamics (WSINDy) and integrates them directly into the framework's PDE assembly pipeline.

---

### Struct Fields & Method Implementations

| Field / Method | Type | Function / Purpose |
| --- | --- | --- |
| `km_expr` | `Symbolics.Num` | Discovered symbolic equation for eddy viscosity of momentum $K_m(\mathbf{z})$. |
| `kh_expr` | `Symbolics.Num` | Discovered symbolic equation for eddy diffusivity of heat $K_h(\mathbf{z})$. |
| `flux_expr` | `Symbolics.Num` | Discovered symbolic equation for boundary layer surface flux $\tau_s$ or $H_s$. |
| `eddy_momentum(c, state)` | Method | Returns symbolic expression `c.km_expr`. |
| `eddy_heat(c, state)` | Method | Returns symbolic expression `c.kh_expr`. |
| `surface_flux(c, state)` | Method | Returns symbolic expression `c.flux_expr`. |

---

### Architectural Role in the Pipeline

```
  Field Campaign Data / Observations
                 │
                 ▼
       WSINDy Discovery Engine
  (Weak-form STRidge + Constrained QP)
                 │
                 ▼
        Discovered Symbolics
   (km_expr, kh_expr, flux_expr)
                 │
                 ▼
          WSINDyClosure
                 │
                 ▼
   ModelingToolkit.jl (build_pde_system)
 (Prognostic Continuous-Time Solvers)

```

1. **Symbolic Integration:** By storing `Symbolics.Num` types rather than hardcoded Julia functions, `WSINDyClosure` allows discovered equations to be compiled into `ModelingToolkit.jl` systems (`build_pde_system`), automatically yielding exact analytical Jacobians for implicit time integrators.
2. **Seamless Interface Compliance:** By overloading `eddy_momentum`, `eddy_heat`, and `surface_flux`, `WSINDyClosure` can be swapped transparently with empirical baselines (such as `MOSTClosure`) or physical scaling laws (`PhysicalSimilarityClosure`) without altering the PDE discretization backends (`MethodOfLinesFD` or `SpectralBLGalerkin`).