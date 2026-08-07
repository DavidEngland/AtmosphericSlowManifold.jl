`src/System/PrognosticPDE.jl` builds a 1D vertical boundary layer PDE system using `ModelingToolkit.jl`. It maps intrinsic manifold state variables ($\eta_1, \eta_2, \eta_3, r, \omega, \chi, \pi_g$) directly into prognostic PDE fields $(u(t,z), v(t,z), \theta(t,z))$, enabling analytical symbolic derivation of closure diffusivities ($K_m, K_h$) and surface fluxes.

---

### Key Mathematical & Physical Components

1. **Prognostic PDE Governing Equations:**

$$\frac{\partial u}{\partial t} = f (v - v_g) + \frac{\partial}{\partial z}\left( K_m \frac{\partial u}{\partial z} \right)$$


$$\frac{\partial v}{\partial t} = -f u + \frac{\partial}{\partial z}\left( K_m \frac{\partial v}{\partial z} \right)$$


$$\frac{\partial \theta}{\partial t} = \frac{\partial}{\partial z}\left( K_h \frac{\partial \theta}{\partial z} \right) + R_Q$$


2. **Manifold Mapping (`manifold_subs`):**
* Boundary layer state variables map continuously to continuous fields: $u(t,z)$, $v(t,z)$, $\theta(t,z)$.
* Shear invariants $\pi_g = \vert{}\partial_z \mathbf{u}\vert{}^2 = (\partial_z u)^2 + (\partial_z v)^2$ and lapse rate $\chi = \partial_z \theta$ are symbolically expanded via `ModelingToolkit.expand_derivatives`.



---

### Upgraded `PrognosticPDE.jl` Implementation

Integrating `default_surface_heat_flux` from `SurfaceBoundary.jl` replaces the static `0.0` lower boundary condition for $\theta$ with dynamic surface heat exchange:

```julia
# src/System/PrognosticPDE.jl
"""
    build_pde_system(closure::AbstractClosure; kwargs...)

Construct a 1D vertical atmospheric boundary layer `ModelingToolkit.PDESystem` for
prognostic momentum (u, v) and potential temperature (theta) transport.
"""
function build_pde_system(
    closure::AbstractClosure;
    z_top::Float64 = 3000.0,
    t_end::Float64 = 86400.0,
    coriolis::Float64 = 1e-4,
    v_geostrophic::Float64 = 8.0,
    radiation::Float64 = 0.0,
    z0_value::Float64 = 0.1,
    u_star_value::Float64 = 0.3,
)
    ModelingToolkit.@parameters t z
    ModelingToolkit.@variables u(..) v(..) theta(..)

    Dt = ModelingToolkit.Differential(t)
    Dz = ModelingToolkit.Differential(z)

    state = ManifoldState()

    # Bind intrinsic manifold symbols to prognostic/diagnostic fields so
    # closures discovered in manifold coordinates are executable in PDE space.
    manifold_subs = Dict(
        state.eta1 => u(t, z),
        state.eta2 => v(t, z),
        state.eta3 => theta(t, z),
        state.r => sqrt(u(t, z)^2 + v(t, z)^2 + 1e-12),
        state.omega => atan(v(t, z), u(t, z) + 1e-12),
        state.chi => Dz(theta(t, z)),
        state.pi_g => Dz(u(t, z))^2 + Dz(v(t, z))^2,
        state.lambdamin => one(z),
        state.u => u(t, z),
        state.v => v(t, z),
        state.theta => theta(t, z),
        state.q => 0.0,
        state.u_star => u_star_value,
        state.z => z,
        state.z0 => z0_value,
    )

    km = ModelingToolkit.expand_derivatives(Symbolics.substitute(eddy_momentum(closure, state), manifold_subs))
    kh = ModelingToolkit.expand_derivatives(Symbolics.substitute(eddy_heat(closure, state), manifold_subs))

    # Surface fluxes evaluated via SurfaceBoundary interface
    tau_s = ModelingToolkit.expand_derivatives(Symbolics.substitute(default_surface_flux(closure, state), manifold_subs))
    q_s = ModelingToolkit.expand_derivatives(Symbolics.substitute(default_surface_heat_flux(closure, state), manifold_subs))

    eqs = [
        Dt(u(t, z)) ~ coriolis * (v(t, z) - v_geostrophic) - Dz(-km * Dz(u(t, z))),
        Dt(v(t, z)) ~ -coriolis * u(t, z) - Dz(-km * Dz(v(t, z))),
        Dt(theta(t, z)) ~ -Dz(-kh * Dz(theta(t, z))) + radiation,
    ]

    bcs = [
        Dz(u(t, 0.0)) ~ tau_s,
        Dz(v(t, 0.0)) ~ 0.0,
        Dz(theta(t, 0.0)) ~ q_s,
        u(t, z_top) ~ v_geostrophic,
        v(t, z_top) ~ 0.0,
        theta(t, z_top) ~ 0.0,
        u(0.0, z) ~ 0.0,
        v(0.0, z) ~ 0.0,
        theta(0.0, z) ~ 0.0,
    ]

    domains = [
        t in Interval(0.0, t_end),
        z in Interval(0.0, z_top),
    ]

    ModelingToolkit.@named pde = ModelingToolkit.PDESystem(
        eqs,
        bcs,
        domains,
        [t, z],
        [u(t, z), v(t, z), theta(t, z)]
    )
    return pde
end

```

---

### Pipeline Integration & Downstream Features

| Component | Functionality |
| --- | --- |
| **MethodOfLines.jl Spatial Discretization** | Converts `PDESystem` into a spatial ODE system ($\partial_t \mathbf{y} = \mathbf{f}(\mathbf{y})$) using finite differences on regular or staggered $z$-grids. |
| **Symbolic Differentiation** | `Symbolics.substitute` preserves $C^\infty$ smooth functions from `SmoothOperators.jl`, ensuring MTK auto-differentiates explicit analytical Jacobians $\mathbf{J} = \frac{\partial \mathbf{f}}{\partial \mathbf{y}}$. |
| **Bifurcation & Continuation** | Discretized ODE system outputs interface directly with `BifurcationKit.jl` / `Continuation.jl` for fast-slow manifold fold tracking. |

---

The `build_pde_system` constructor now fully integrates dynamic lower surface boundary conditions for both momentum ($\tau_s$) and potential temperature ($q_s$) within the `ModelingToolkit.jl` symbolic pipeline.

---

### Updated Boundary Condition Matrix

| Field | Spatial Boundary | Boundary Form ($z = 0$) | Boundary Form ($z = z_{\text{top}}$) | Physical Mechanism |
| --- | --- | --- | --- | --- |
| **$u(t, z)$** | Lower / Upper | $\frac{\partial u}{\partial z} = \tau_s$ | $u(t, z_{\text{top}}) = v_g$ | Surface stress drive / Geostrophic velocity match |
| **$v(t, z)$** | Lower / Upper | $\frac{\partial v}{\partial z} = 0$ | $v(t, z_{\text{top}}) = 0$ | Zero cross-wind surface shear / Free atmosphere relaxation |
| **$\theta(t, z)$** | Lower / Upper | $\frac{\partial \theta}{\partial z} = q_s$ | $\theta(t, z_{\text{top}}) = 0$ | Surface heat exchange ($q_s$) / Upper thermal boundary |

---

### System Implications

1. **Stratification & Stability Coupling:**
With dynamic lower boundary conditions $\frac{\partial \theta}{\partial z} = q_s$, surface heat fluxes directly modulate vertical thermal gradients $\chi = \frac{\partial \theta}{\partial z}$ in the boundary layer. This enables realistic diurnal heating/cooling simulations where Monin–Obukhov stability parameters ($\zeta = z/L$) evolve dynamically over time.
2. **Symbolic Differentiation & Jacobian Construction:**
Because $q_s$ is ingested via `default_surface_heat_flux` and expanded through `Symbolics.substitute`, `ModelingToolkit.expand_derivatives` preserves $C^\infty$ smoothness. When discretized using `MethodOfLines.jl`, exact analytical Jacobians ($\mathbf{J} = \frac{\partial \mathbf{f}}{\partial \mathbf{y}}$) remain available for stiff ODE solvers and fast-slow manifold continuation.
3. **Zero-Overhead Numeric Fallback:**
In non-thermal or unheated test cases where `q_s` evaluates to `0.0`, the system automatically defaults to zero lower heat flux without requiring separate PDE setup pathways.

---

`build_pde_system` constructs a `ModelingToolkit.PDESystem` representing 1D vertical atmospheric boundary layer dynamics for momentum ($u, v$) and thermal transport ($\theta$), fully integrated with `ManifoldState` coordinate transformations and `SurfaceBoundary` flux dispatch.

---

### System Architecture & Field Mapping

| Component | Symbolic Representation | Operational Mapping |
| --- | --- | --- |
| **Prognostic Variables** | `u(t, z)`, `v(t, z)`, `theta(t, z)` | Primary velocity components ($u, v$) and potential temperature ($\theta$). |
| **Shear Invariant ($\pi_g$)** | `Dz(u(t, z))^2 + Dz(v(t, z))^2` | Local squared vertical velocity shear $\vert{}\partial_z \mathbf{u}\vert{}^2$. |
| **Lapse Rate ($\chi$)** | `Dz(theta(t, z))` | Local vertical potential temperature gradient $\partial_z \theta$. |
| **Speed & Direction ($r, \omega$)** | `sqrt(u^2 + v^2 + 1e-12)`, `atan(v, u + 1e-12)` | Regulated horizontal wind velocity magnitude and polar flow angle. |
| **Eddy Diffusivities** | `km`, `kh` | Expanded analytical expressions for momentum ($K_m$) and heat ($K_h$) transport derived from the active `closure`. |

---

### Governing PDE Mechanics & Boundary Specifications

#### 1. Conservation Laws

$$\frac{\partial u}{\partial t} = f (v - v_g) + \frac{\partial}{\partial z}\left( K_m \frac{\partial u}{\partial z} \right)$$

$$\frac{\partial v}{\partial t} = -f u + \frac{\partial}{\partial z}\left( K_m \frac{\partial v}{\partial z} \right)$$

$$\frac{\partial \theta}{\partial t} = \frac{\partial}{\partial z}\left( K_h \frac{\partial \theta}{\partial z} \right) + R_Q$$

#### 2. Surface & Domain Conditions ($z \in [0, z_{\text{top}}]$)

* **Surface Layer ($z = 0$ Neumann):**

$$\left. \frac{\partial u}{\partial z} \right\vert{}_{z=0} = \tau_s, \qquad \left. \frac{\partial v}{\partial z} \right\vert{}_{z=0} = 0, \qquad \left. \frac{\partial \theta}{\partial z} \right\vert{}_{z=0} = q_s$$


* **Free Atmosphere ($z = z_{\text{top}}$ Dirichlet):**

$$u(t, z_{\text{top}}) = v_g, \qquad v(t, z_{\text{top}}) = 0, \qquad \theta(t, z_{\text{top}}) = 0$$


* **Initial Profile ($t = 0$):**

$$u(0, z) = 0, \qquad v(0, z) = 0, \qquad \theta(0, z) = 0$$



---

### Numerical Solver Compatibility

* **`MethodOfLines.jl` Ready:** The symbolic output `PDESystem` can be directly spatial-discretized using `MOLify(grid_spec)` into an `ODEProblem` or `DAEProblem`.
* **Exact Analytical Jacobians:** Because `Symbolics.substitute` and `ModelingToolkit.expand_derivatives` preserve smooth functions from $C^\infty$ closures, automatic symbolic differentiation provides analytical system Jacobians $\mathbf{J} = \frac{\partial \mathbf{f}}{\partial \mathbf{y}}$ for stiff implicit time integrators (such as `TRBDF2` or `KenCarp4`).