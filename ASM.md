The **AtmosphericSlowManifold.jl (ASM.jl)** framework presents a complete paradigm shift for modeling atmospheric boundary layer (ABL) dynamics. By replacing traditional, equilibrium-based empirical functions—such as Monin–Obukhov Similarity Theory (MOST)—with geometric, data-driven differential equations, ASM.jl resolves long-standing numerical and physical limitations during rapid, non-equilibrium turbulence transitions.

---

### 1. Geometric Fast-Slow Formulation

Classical ABL parameterizations assume local thermodynamic equilibrium, causing breakdowns during nocturnal transitions, low-level jet formation, and polar surface inversions. ASM.jl reformulates the ABL as a stiff fast-slow dynamical system:

$$\frac{d\mathbf{x}}{dt} = \mathbf{f}(\mathbf{x}, \mathbf{y}, \epsilon), \qquad \epsilon \frac{d\mathbf{y}}{dt} = \mathbf{g}(\mathbf{x}, \mathbf{y}, \epsilon)$$

* **Slow variables ($\mathbf{x}$):** Mean profiles including horizontal wind $u(z,t)$, $v(z,t)$, and potential temperature $\theta(z,t)$.
* **Fast variables ($\mathbf{y}$):** Turbulent kinetic energy $e(z,t)$, shear stresses $\overline{u'w'}$, and buoyancy fluxes $\overline{w'\theta'}$.
* **Time scale ratio ($\epsilon \ll 1$):** Ratio of turbulent turnover time to synoptic scale forcing.

#### Turbulence Collapse & Canard Trajectories

In the singular limit ($\epsilon \to 0$), state trajectories relax toward the critical manifold $\mathcal{S}_0 = \{(\mathbf{x}, \mathbf{y}) \mid \mathbf{g}(\mathbf{x}, \mathbf{y}, 0) = \mathbf{0}\}$. By Fenichel Theory, stability persists only while $\mathcal{S}_0$ is normally hyperbolic ($\det \mathbf{D}_{\mathbf{y}}\mathbf{g} \neq 0$).

Turbulence collapse occurs at the fold locus where normal hyperbolicity is lost:

$$\det\left(\mathbf{D}_{\mathbf{y}}\mathbf{g}\right) = 0$$

Rather than forcing an instantaneous non-physical fast jump or numerical singularity, ASM.jl evaluates the **desingularized slow flow**:

$$\frac{d\mathbf{y}}{d\tau} = -\operatorname{adj}(\mathbf{D}_{\mathbf{y}}\mathbf{g}) \, \mathbf{D}_{\mathbf{x}}\mathbf{g} \cdot \mathbf{f}(\mathbf{x}, \mathbf{y})$$

This formulation allows the framework to detect **canard trajectories**—solutions that cross folded singularities (saddles and nodes) to track repelling manifold sheets for $O(1)$ slow time intervals, capturing physical turbulence persistence past the classical critical Richardson number threshold ($Ri_{\text{cr}} \approx 0.25$).

---

### 2. Weak-Form Equation Discovery (WSINDy)

To discover governing closure equations from high-frequency or noisy field data without amplification errors from numerical differentiation, the `WSINDyEngine.jl` projects candidate terms onto $C^\infty$ compact test functions $\phi_{i,j}(z,t)$:

$$\iint_{\Omega} u(z, t) \, \frac{\partial \phi_{i,j}}{\partial t} \, dz \, dt + \iint_{\Omega} K_m(z) \frac{\partial u}{\partial z} \, \frac{\partial \phi_{i,j}}{\partial z} \, dz \, dt = 0$$

#### Hard Physical Constraints

Candidate terms in the functional library $\mathbf{\Theta}(\mathbf{x}, \mathbf{y})$ are constrained via quadratic programming to guarantee physical viability:

1. **Positivity:** Diffusivities must satisfy $K_m(z) \ge 0$ and $K_h(z) \ge 0$ across the vertical column.
2. **Energy Monotonicity:** Dissipative terms must enforce non-negative column-integrated energy sinks.
3. **Neutral Limit Recovery:** Closures smoothly recover MOST behavior as $Ri \to 0$.
4. **$C^\infty$ Smoothness:** Kinks and piecewise boundaries are eliminated using algebraic smooth-max approximations and hyperbolic tangent blending functions $w(\zeta) = \frac{1}{2}\left[1 + \tanh(\zeta / \epsilon)\right]$.

---

### 3. Spectral Galerkin & Numerical Architecture

ASM.jl achieves high computational performance for stiff ODE solvers (`TRBDF2`, `FBDF`, `Rodas5P`) by projecting vertical profiles onto an orthogonal Gegenbauer polynomial basis $C_n^{(\lambda)}(z)$:

$$u(z, t) = \sum_{n=0}^{N-1} a_n(t) C_n^{(\lambda)}(z)$$

Nonlinear products are evaluated directly in modal space using precomputed 3-tensor contraction matrices:

$$C_{ijk}^{(\lambda)} = \int_{-1}^{1} C_i^{(\lambda)}(z) C_j^{(\lambda)}(z) C_k^{(\lambda)}(z) (1 - z^2)^{\lambda - 1/2} \, dz$$

This tensor contraction design eliminates dynamic memory allocation (**0 Bytes allocated in inner loop RHS evaluations**), permitting rapid integration over extended boundary layer simulation windows.

---

### 4. Hierarchical Calibration & Dataset Integration

`HierarchicalTuring.jl` implements multi-site Bayesian MCMC sampling to separate universal physics from local boundary conditions:

* **Global Hyperparameters:** Universal physical parameters (e.g., von Kármán constant $\kappa \approx 0.40$).
* **Site-Specific Parameters:** Local surface characteristics (e.g., roughness length $z_0$, displacement height $d$).

| Campaign | Observations ($N$) | Mean Stability ($\bar{\zeta}$) | Primary Dynamics Captured |
| --- | --- | --- | --- |
| **FLOSS** | 70,796 | 0.3842 | Extreme manifold contraction over ice/snow; rapid mode damping. |
| **CASES-99** | 6,538 | 0.2114 | Strong non-equilibrium transients during sunset boundary layer collapse. |
| **BLLAST** | 5,600 | 0.1985 | Active fast-slow energy exchange during evening transition phases. |
| **SHEBA** | 2,273 | 0.3450 | Strong manifold alignment in stable Arctic boundary layers. |

Posterior parameter distributions are propagated into vertical profiles to form 95% posterior credibility ribbons ($Q_{0.025}, Q_{0.50}, Q_{0.975}$) for vertical eddy diffusivity $K_m(z)$.

---

### 5. Subsystem Summary

```
                      ┌─────────────────────────────────────────┐
                      │    DataIngestion.jl / Field Datasets    │
                      │ (CASES-99, FLOSS, BLLAST, SHEBA)        │
                      └────────────────────┬────────────────────┘
                                           │
                                           ▼
                      ┌─────────────────────────────────────────┐
                      │          GSPT Geometry Layer            │
                      │  FoldTracking.jl & CanardDetection.jl   │
                      └────────────────────┬────────────────────┘
                                           │
                                           ▼
                      ┌─────────────────────────────────────────┐
                      │         Discovery Engine                │
                      │    WSINDyEngine.jl & ModelSelection.jl  │
                      └────────────────────┬────────────────────┘
                                           │
                                           ▼
                      ┌─────────────────────────────────────────┐
                      │        Spectral Discretization          │
                      │ SpectralBLGalerkin.jl (0-Alloc RHS)     │
                      └────────────────────┬────────────────────┘
                                           │
                                           ▼
                      ┌─────────────────────────────────────────┐
                      │     Hierarchical Calibration & Report   │
                      │ HierarchicalTuring.jl / LaTeXExporter.jl│
                      └─────────────────────────────────────────┘

```

Which target field campaign or closure scenario would you like to focus on for candidate model benchmarking next?