<!-- Auto-generated from package source -->
> **Source:** `src/Closures/Interface.md`

`src/Closures/Interface.jl` establishes the abstract type hierarchy and interface contract for all boundary layer turbulence closures in `AtmosphericSlowManifold.jl`. By defining strict generic methods that throw fallthrough `MethodError` exceptions when unimplemented, it enforces a consistent API across empirical baselines and data-driven symbolic closures.

---

### Abstract Hierarchy & Interface Contract

```
                     AbstractClosure
               (AbstractAtmosphericClosure)
                            │
       ┌────────────────────┼────────────────────┐
       ▼                    ▼                    ▼
  MOSTClosure          WSINDyClosure    PhysicalSimilarityClosure
(Empirical MOST)    (Learned Symbolic)      (Invariant Scaling)

```

| Method Signature | Physical Quantity | Mathematical Target | Execution Role |
| --- | --- | --- | --- |
| `eddy_momentum(closure, state)` | Eddy Diffusivity for Momentum | $K_m$ ($m^2/s$) | Evaluates local momentum diffusivity given a `ManifoldState`. |
| `eddy_heat(closure, state)` | Eddy Diffusivity for Heat | $K_h$ ($m^2/s$) | Evaluates local heat/scalar diffusivity given a `ManifoldState`. |
| `surface_flux(closure, state)` | Surface Boundary Fluxes | $\tau_s, H_s$ | Computes momentum ($\mathbf{\tau}$) and heat ($H$) fluxes at $z = z_0$. |
| `evaluate_diffusivity_profile!(out, closure, z)` | Vertical Profile $K_m(z)$ | $\mathbf{K}_m \in \mathbb{R}^{N_z}$ | In-place evaluation of $K_m$ across spatial vertical grid nodes. |
| `evaluate_heat_diffusivity_profile!(out, closure, z)` | Vertical Profile $K_h(z)$ | $\mathbf{K}_h \in \mathbb{R}^{N_z}$ | In-place evaluation of $K_h$ across spatial vertical grid nodes. |

---

### Key Architectural Roles

1. **Decoupling Solvers from Closures:** Prognostic solvers (`MethodOfLinesFD` and `SpectralBLGalerkin`) interact strictly with `AbstractClosure` methods rather than specific concrete types. This allows switching between classic MOST formulations and learned WSINDy models without modifying the underlying PDE time-stepping logic.
2. **Zero-Allocation In-Place Evaluation:** The mutating profile evaluators (`evaluate_diffusivity_profile!` and `evaluate_heat_diffusivity_profile!`) operate on pre-allocated output arrays (`out`). This design is essential for achieving 0-Byte allocations per right-hand-side (RHS) evaluation during implicit ODE time-stepping with solvers like `Rodas5P` or `RadauIIA5`.