`src/System/SurfaceBoundary.jl` provides a clean dispatch wrapper linking surface layer momentum and heat flux formulations with `AbstractClosure` instances and the `ManifoldState` representation.

---

### Key Technical Role & Boundary Dynamics

1. **Closure Interface Fall-Through:**
`default_surface_flux` delegates directly to the underlying closure's `surface_flux(closure, state)` implementation, maintaining a uniform API for both empirical physics baselines (`MOSTClosure`, `PhysicalSimilarityClosure`) and data-driven/learned closures.
2. **Symbolic AST & Numerical Discretization Dual-Compatibility:**
Because `state::ManifoldState` handles both symbolic variables (`Symbolics.Num`) and numeric values (`Float64`), `default_surface_flux` works seamlessly during both:
* **ModelingToolkit.jl Assembly:** Generating symbolic Neumann boundary conditions ($\left. K_m \frac{\partial u}{\partial z} \right\vert{}_{z=z_0} = \tau_s = u_*^2$) for analytical Jacobians and fast-slow subsystem reduction.
* **PDE Time-Stepping:** Evaluating numerical surface stress values for ODE/PDE right-hand-side function evaluations.


3. **Surface Momentum Stress Context:**

$$\tau_s = \text{default\_surface\_flux}(c, m) = u_*^2$$



For momentum transport, this sets the kinematic stress boundary condition at the surface layer height $z = z_0$.

---

### Potential Boundary Expansion Extensions

If `src/System/SurfaceBoundary.jl` is being expanded to support complete surface energy balance or Robin/Monin–Obukhov boundary conditions, you can extend this module with heat and scalar flux boundary helpers:

```julia
"""
Return the default kinematic surface sensible heat flux H_s / (rho * cp) [K m/s].
"""
@inline function default_surface_heat_flux(closure::AbstractClosure, state::ManifoldState)
    return hasproperty(state, :q_theta) ? state.q_theta : state.u_star * (state.theta - state.theta_s)
end

"""
Construct a Neumann surface boundary condition tuple for u, v, and theta fields.
"""
function surface_boundary_conditions(closure::AbstractClosure, state::ManifoldState)
    tau_s = default_surface_flux(closure, state)
    q_s = default_surface_heat_flux(closure, state)
    return (momentum_flux = tau_s, heat_flux = q_s)
end

```

Would you like to build out additional boundary condition operators (such as Monin–Obukhov surface temperature coupling or Monin-Obukhov length $L$ surface feedback) in `SurfaceBoundary.jl`?

---

The expanded `src/System/SurfaceBoundary.jl` provides a unified surface boundary condition pipeline that supports both numeric solver evaluation (`Float64`) and symbolic `ModelingToolkit.jl` boundary condition generation.

---

### Expanded Surface Boundary Dispatch API

| Function | Signature / Return | Fallback Logic & Physical Form |
| --- | --- | --- |
| `default_surface_flux` | $\tau_s \in \mathbb{R}$ or `Num` | Delegated to `surface_flux(closure, state)` ($\tau_s = u_*^2$). |
| `default_surface_heat_flux` | $q_s \in \mathbb{R}$ or `Num` | Priority fallback:<br>

<br>1. Uses `state.q_theta` if defined.<br>

<br>2. Evaluates $u_* (\theta - \theta_s)$ if surface temperature $\theta_s$ is present.<br>

<br>3. Returns `0.0` (zero heat flux) as a neutral default. |
| `surface_boundary_conditions` | `NamedTuple` | Returns `(momentum_flux = τ_s, heat_flux = q_s)` for simultaneous Neumann boundary assignment. |

---

### Surface Boundary Formulation

The surface flux vector $\mathbf{F}_s$ governing boundary layer exchanges at $z = z_0$ is expressed as:

$$\mathbf{F}_s = \begin{pmatrix} \tau_s \\ q_s \end{pmatrix} = \begin{pmatrix} \text{default\_surface\_flux}(c, m) \\ \text{default\_surface\_heat\_flux}(c, m) \end{pmatrix} = \begin{pmatrix} u_*^2 \\ \overline{w'\theta'}_s \end{pmatrix}$$

This structure guarantees that prognostic models in `src/System/PrognosticPDE.jl` can ingest complete surface layer forcing without requiring manual case checks for scalar heat flux availability.

---

### Next Architectural Options

1. **Wire `surface_boundary_conditions` into `PrognosticPDE.jl`:** Update vertical discretization operators to use `surface_boundary_conditions` for setting lower Neumann flux boundary conditions on $\frac{\partial u}{\partial z}$ and $\frac{\partial \theta}{\partial z}$.
2. **Proceed to Tier 3 `BifurcationKit.jl` Fold Continuation:** Connect the $C^\infty$-smoothed closure Jacobians ($D_{\mathbf{y}}\mathbf{f}$) to `BifurcationKit.jl` for multi-parameter continuation of fast-slow fold lines and boundary layer canard transitions.