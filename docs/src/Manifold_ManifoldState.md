<!-- Auto-generated from package source -->
> **Source:** `src/Manifold/ManifoldState.md`

`src/Manifold/ManifoldState.jl` defines the core symbolic data structure that bridges physical atmospheric state variables with intrinsic fast-slow manifold coordinates. By declaring all fields as `ModelingToolkit.Num`, it provides a symbolic substrate for constructing manifold mappings, computing GSPT diagnostics, and generating `ModelingToolkit.jl` equations.

### Field Mapping & Physical Significance

| Variable Category | Field Name | Symbolic Variable | Physical / Geometric Interpretation |
| --- | --- | --- | --- |
| **Intrinsic Coordinates** | `r` | $R$ | Modal amplitude representing energetic turbulence scale. |
|  | `omega` | $\Omega$ | Modal phase tracking non-equilibrium dynamics. |
|  | `chi` | $\chi$ | Local manifold curvature. |
|  | `pi_g` | $\Pi_G$ | Conductive ground-coupling parameter. |
|  | `lambdamin` | $\lambda_{\min}$ | Minimum eigenvalue of fast Jacobian $D_{\mathbf{y}}\mathbf{g}$; tracks fold proximity ($\lambda_{\min} \to 0$). |
|  | `eta1`, `eta2`, `eta3` | $\eta_1, \eta_2, \eta_3$ | Reduced slow/manifold coordinate projections. |
| **Physical Observables** | `u`, `v` | $u, v$ | Zonal and meridional wind velocity components ($m/s$). |
|  | `theta` | $\theta$ | Potential temperature ($K$). |
|  | `q` | $q$ | Specific humidity ($kg/kg$). |
| **Boundary & Geometry** | `u_star` | $u_*$ | Friction velocity ($m/s$). |
|  | `z` | $z$ | Height above ground ($m$). |
|  | `z0` | $z_0$ | Surface roughness length ($m$). |

---

### Key Architectural Role

1. **Symbolic Substrate:** The constructor `ManifoldState(; name = :manifold)` invokes `ModelingToolkit.@variables` to generate unbound symbolic variables. This enables downstream modules (`Closures`, `Discovery`, and `PrognosticPDE`) to assemble mathematical expressions symbolically prior to numeric discretization.
2. **Coordinate Decoupling:** Storing both physical observables ($u, v, \theta, q$) and manifold intrinsic variables ($R, \Omega, \chi, \Pi_G$) within a single structure enables direct transformation between physical state space and geometric space during WSINDy feature library evaluation and GSPT fold detection.