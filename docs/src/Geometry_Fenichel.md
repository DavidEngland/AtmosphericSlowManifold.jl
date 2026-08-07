<!-- Auto-generated from package source -->
> **Source:** `src/Geometry/Fenichel.md`

`src/Geometry/Fenichel.jl` implements normal hyperbolicity diagnostics based on **Fenichel’s Theorem**. It quantifies the persistence and stability of slow invariant manifolds under small turbulent perturbations ($\varepsilon > 0$) by analyzing the spectrum of the fast subsystem Jacobian $D_{\mathbf{y}}\mathbf{g}$.

---

### Mathematical Foundation: Fenichel Normal Hyperbolicity

According to Fenichel’s First Theorem, a critical manifold $\mathcal{S}_0 = \{(\mathbf{x}, \mathbf{y}) \mid \mathbf{g}(\mathbf{x}, \mathbf{y}) = \mathbf{0}\}$ persists as a smooth, locally invariant slow manifold $\mathcal{S}_\varepsilon$ for $0 < \varepsilon \ll 1$ **if and only if** $\mathcal{S}_0$ is normally hyperbolic.

Normal hyperbolicity requires that every eigenvalue $\lambda_i$ of the fast subsystem Jacobian $D_{\mathbf{y}}\mathbf{g}$ satisfies:

$$\operatorname{Re}(\lambda_i) \neq 0 \quad \forall \, i \in \{1, \dots, n\}$$

* **Attracting Manifold:** $\operatorname{Re}(\lambda_i) < 0 \; \forall i$ (turbulent fluctuations rapidly decay onto the slow mean boundary layer profile).
* **Loss of Normal Hyperbolicity:** $\min_i \vert{}\operatorname{Re}(\lambda_i)\vert{} \to 0$ (occurs at fold lines, Hopf bifurcations, or points of turbulence collapse).

---

### Structure & Functionality Breakdown

#### 1. Data Structure: `HyperbolicityReport`

Acts as a diagnostic payload containing spectral metrics computed at a single state space point:

| Field | Type | Mathematical Definition | Physical / Numerical Meaning |
| --- | --- | --- | --- |
| `eigvals` | `Vector{ComplexF64}` | $\lambda_i \in \operatorname{spec}(D_{\mathbf{y}}\mathbf{g})$ | Full complex eigenvalue spectrum of the fast subsystem. |
| `min_abs_real` | `Float64` | $\sigma_1 = \min_i \Vert{}\operatorname{Re}(\lambda_i)\Vert{}$ | Rate of attraction/repulsion normal to $\mathcal{S}_0$. Tracks fold proximity ($\sigma_1 \to 0$). |
| `spectral_gap` | `Float64` | $\Delta = \sigma_2 - \sigma_1$ | Separation between dominant and sub-dominant fast dynamics modes. |
| `epsilon0` | `Float64` | $\varepsilon_0 \propto \sigma_1$ | Upper bound on perturbation size $\varepsilon$ for manifold persistence. |
| `is_hyperbolic` | `Bool` | $\sigma_1 > \text{gap\_tol}$ | Boolean flag indicating whether the slow manifold remains persistent. |

---

#### 2. Multiple Dispatch Engine: `fenichel_metrics`

The diagnostic engine provides three overloaded interfaces to process different representation layers:

```
                         Input Types
                             │
     ┌───────────────────────┼───────────────────────┐
     ▼                       ▼                       ▼
Explicit Matrix        Black-Box Vector Field    Compiled Model
J::AbstractMatrix      f_fast(x, y, p)          model::JacobianModel
     │                       │                       │
     │                       ▼                       │
     │            Central Finite Differences          │
     │         finite_difference_jacobian_y()        │
     │                       │                       │
     └───────────────────────┼───────────────────────┘
                             │
                             ▼
                    eigenvalues spec(J)
                             │
                             ▼
                   HyperbolicityReport

```

* **Direct Matrix:** `fenichel_metrics(J; gap_tol)` extracts eigenvalues directly from pre-computed numerical matrices.
* **Finite-Difference Fallback:** `fenichel_metrics(f_fast, x, y; params, gap_tol)` uses 2nd-order central differences to approximate $D_{\mathbf{y}}\mathbf{g}$ when symbolic definitions are missing.
* **Symbolic Compiled Kernel:** `fenichel_metrics(model::JacobianModel, u_fast, u_slow, p; gap_tol)` uses pre-compiled C-speed function evaluations (`evaluate_jacobian`) for high-throughput time-stepping and trajectory diagnostics.

---

#### 3. Surface Spatial Mapping: `hyperbolicity_profile`

```julia
function hyperbolicity_profile(surface::CriticalManifoldSurface, f_fast; ...)

```

Maps `fenichel_metrics` across an entire `CriticalManifoldSurface`. It generates a vector of `HyperbolicityReport` objects along a vertical boundary layer profile or spatial path, isolating regions where static stability forces $D_{\mathbf{y}}\mathbf{g}$ toward singular regimes ($\sigma_1 \to 0$).

---

### Physical Utility in Boundary Layer Modeling

1. **Quantifying Turbulence Breakdown:** Provides an explicit bound on when standard turbulence closure approximations fail due to the destruction of normal attraction.
2. **Step Size Control:** Guides adaptive ODE/PDE integrators by lowering step sizes near non-hyperbolic regions ($\sigma_1 \le 10^{-6}$) to properly resolve fast-slow cross-talk and canard transitions.