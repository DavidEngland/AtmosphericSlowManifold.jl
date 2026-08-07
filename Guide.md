# AtmosphericSlowManifold.jl: A Gentle Guide to Data-Driven Atmospheric Dynamics & Turbulence Closure

`AtmosphericSlowManifold.jl` (ASM.jl) is a high-performance Julia framework designed to model, analyze, and discover atmospheric boundary layer (ABL) dynamics. By unifying **Geometric Singular Perturbation Theory (GSPT)** with **Weak-form Sparse Identification of Non-linear Dynamics (WSINDy)**, ASM.jl replaces empirical, steady-state turbulence closures with data-driven, physically constrained differential equations.

---

## 1. Core Philosophy & Architectural Overview

Classical turbulence closure models—such as Monin–Obukhov Similarity Theory (MOST)—rely on empirical stability functions $\phi_m(\zeta)$ and $\phi_h(\zeta)$ that assume local thermodynamic equilibrium. These formulations breakdown during rapid non-equilibrium transitions, such as sunset boundary layer collapse, low-level jet formation, and extreme stability over polar ice sheets.

ASM.jl addresses these failures by decoupling observation space, intrinsic manifold geometry, and spatial discretization.

```
                  ┌─────────────────────────────────────────┐
                  │   Raw Field Observations & Tower CSVs   │
                  │   (CASES-99, FLOSS, BLLAST, SHEBA)      │
                  └────────────────────┬────────────────────┘
                                       │
                                       ▼
                  ┌─────────────────────────────────────────┐
                  │    Data Ingestion & Gegenbauer Bases    │
                  │        (src/Observation/)               │
                  └────────────────────┬────────────────────┘
                                       │
            ┌──────────────────────────┴──────────────────────────┐
            ▼                                                     ▼
┌───────────────────────┐                             ┌───────────────────────┐
│ Fast-Slow GSPT Engine │                             │  Weak-Form WSINDy     │
│ (src/Geometry/)       │                             │  (src/Discovery/)     │
│ • Critical Manifolds  │                             │ • Space-time test fn  │
│ • Fold Locus det=0    │                             │ • Hard Positivity QP  │
│ • Canard Tracking     │                             │ • Pareto AIC/BIC      │
└───────────┬───────────┘                             └───────────┬───────────┘
            │                                                     │
            └──────────────────────────┬──────────────────────────┘
                                       │
                                       ▼
                  ┌─────────────────────────────────────────┐
                  │      C^∞ Differentiable Closures        │
                  │        (src/Closures/)                  │
                  └────────────────────┬────────────────────┘
                                       │
            ┌──────────────────────────┴──────────────────────────┐
            ▼                                                     ▼
┌───────────────────────┐                             ┌───────────────────────┐
│ Discretization Engine │                             │ Bayesian Calibration  │
│ (src/Discretization/) │                             │ (src/Calibration/)    │
│ • Stretched FD        │                             │ • MLE / VI / MCMC     │
│ • 3-Tensor Spectral   │                             │ • Hierarchical Turing │
│   Galerkin (0 Alloc)  │                             │ • 95% K_m(z) Ribbons  │
└───────────┬───────────┘                             └───────────┬───────────┘
            │                                                     │
            └──────────────────────────┬──────────────────────────┘
                                       │
                                       ▼
                  ┌─────────────────────────────────────────┐
                  │   Automated LaTeX & Report Exporters    │
                  │   (reports/ & LaTeXExporter.jl)         │
                  └─────────────────────────────────────────┘

```

---

## 2. Framework Subsystems & Module Map

The architecture of ASM.jl is partitioned into seven primary functional domains, located under `src/`:

| Layer | Primary Submodule | Key Types & Functions | Scientific Role |
| --- | --- | --- | --- |
| **Observation** | `DataIngestion.jl` | `ObservationTable`, `read_tower_csv`, `project_to_gegenbauer` | Ingests CSV/NetCDF field campaign data, resolves variable aliases ($H, u_*$), and projects profiles onto Gegenbauer polynomials. |
| **Geometry** | `FoldTracking.jl`<br>

<br>`CanardDetection.jl` | `ManifoldState`, `FoldConstraint`, `compute_jacobian` | Tracks critical slow manifolds $\mathcal{S}_0$, calculates fast subsystem Jacobians $\mathbf{D}_{\mathbf{y}}\mathbf{g}$, detects fold points, and measures normal hyperbolicity. |
| **Discovery** | `WSINDyEngine.jl`<br>

<br>`ModelSelection.jl` | `build_weak_library`, `fit_wsindy_jump`, `compute_pareto_front` | Constructs noise-robust weak-form integrations, enforces physical linear inequality bounds via JuMP.jl, and evaluates AIC/BIC Pareto trade-off curves. |
| **Closures** | `MOSTClosure.jl`<br>

<br>`WSINDyClosure.jl` | `AbstractAtmosphericClosure`, `MOSTClosure`, `WSINDyClosure` | Houses $C^\infty$-differentiable eddy diffusivity profiles ($K_m, K_h$). Eliminates step-discontinuities using hyperbolic blending functions. |
| **Discretization** | `StretchedGrid.jl`<br>

<br>`SpectralBLGalerkin.jl` | `MethodOfLinesFD`, `SpectralBLGalerkin`, `generate_stretched_grid` | Offers stretched 1D vertical finite differences and zero-allocation Gegenbauer-Galerkin spectral expansions using 3-tensor contractions ($C_{ijk}^{(\lambda)}$). |
| **Calibration** | `HierarchicalTuring.jl`<br>

<br>`UncertaintyPropagation.jl` | `MaximumLikelihood`, `BayesianMCMC`, `evaluate_profile_uncertainty` | Performs deterministic OLS/Ridge estimation, Turing.jl MCMC sampling, multi-site hierarchical pooling, and 95% spatial profile credibility ribbon evaluations. |
| **Reporting** | `LaTeXExporter.jl`<br>

<br>`ExportUtilities.jl` | `to_latex`, `latex_term_table`, `latex_site_summary_table` | Converts discovered symbolic ASTs into publication-ready LaTeX differential equations and tabular environments. |

---

## 3. Mathematical Foundations

### 3.1 Geometric Singular Perturbation Theory (GSPT)

The atmospheric boundary layer is modeled as a fast-slow dynamical system:

$$\frac{d\mathbf{x}}{dt} = \mathbf{f}(\mathbf{x}, \mathbf{y}, \epsilon), \qquad \epsilon \frac{d\mathbf{y}}{dt} = \mathbf{g}(\mathbf{x}, \mathbf{y}, \epsilon)$$

* $\mathbf{x} \in \mathbb{R}^m$: Slow mean variables (wind speed $u$, potential temperature $\theta$).
* $\mathbf{y} \in \mathbb{R}^n$: Fast turbulent variables (turbulent kinetic energy $e$, Reynolds stresses $\overline{u'w'}$).
* $\epsilon \ll 1$: Ratio of turbulent timescale to synoptic evolution timescale.

#### Critical Manifold & Normal Hyperbolicity

In the singular limit ($\epsilon \to 0$), fast transients decay onto the **critical manifold**:

$$\mathcal{S}_0 = \left\{ (\mathbf{x}, \mathbf{y}) \in \mathbb{R}^{m+n} \mid \mathbf{g}(\mathbf{x}, \mathbf{y}, 0) = \mathbf{0} \right\}$$

By **Fenichel's Theorem**, $\mathcal{S}_0$ persists as a smooth invariant manifold $\mathcal{S}_\epsilon$ under small perturbations ($\epsilon > 0$) provided it is **normally hyperbolic**—meaning all eigenvalues $\lambda_i$ of the fast Jacobian $\mathbf{D}_{\mathbf{y}}\mathbf{g}$ have non-zero real parts:

$$\operatorname{Re}\left(\lambda_i\left(\mathbf{D}_{\mathbf{y}}\mathbf{g}\right)\right) \neq 0 \quad \forall i \in \{1, \dots, n\}$$

#### Turbulence Collapse & Adjugate Desingularization

Turbulence collapse occurs at the **fold locus**, where normal hyperbolicity is lost:

$$\det\left(\mathbf{D}_{\mathbf{y}}\mathbf{g}(\mathbf{z})\right) = 0$$

To integrate trajectories across singular fold lines, ASM.jl applies **adjugate desingularization** with rescaled time $d\tau = dt / \det(\mathbf{D}_{\mathbf{y}}\mathbf{g})$:

$$\frac{d\mathbf{z}}{d\tau} = \operatorname{adj}\left(\mathbf{D}_{\mathbf{y}}\mathbf{g}\right) \cdot \mathbf{f}(\mathbf{z})$$

This formulation resolves **folded singularities** (saddles, nodes, and foci), allowing the detection of **canard trajectories**—solutions that follow repelling (unstable) manifold sheets without immediate turbulent breakdown.

---

### 3.2 Weak-Form WSINDy Discovery Engine

To avoid amplifying observational noise through direct numerical differentiation, the discovery engine projects governing PDEs against $C^\infty$ test functions $\phi_{i,j}(z, t)$:

$$\iint_{\Omega} \mathcal{L}[u](z, t) \, \phi_{i,j}(z, t) \, dz \, dt = 0$$

Integrating by parts transfers spatial and temporal derivatives directly onto $\phi_{i,j}$:

$$\iint_{\Omega} u(z, t) \, \frac{\partial \phi_{i,j}}{\partial t} \, dz \, dt + \iint_{\Omega} K_m(z) \frac{\partial u}{\partial z} \, \frac{\partial \phi_{i,j}}{\partial z} \, dz \, dt = 0$$

This constructs a linear algebraic system $\mathbf{G}\boldsymbol{\xi} = \mathbf{b}$, solved via sequentially thresholded ridge regression (`STRidge`) or constrained quadratic programming (`ConstrainedQP` using JuMP.jl):

$$\min_{\boldsymbol{\xi}} \Vert{}\mathbf{G}\boldsymbol{\xi} - \mathbf{b}\Vert{}_2^2 + \lambda \Vert{}\boldsymbol{\xi}\Vert{}_1 \quad \text{subject to} \quad \mathbf{C}\boldsymbol{\xi} \ge \mathbf{d}$$

#### Physical Hard Constraints

1. **Positivity:** $K_m(z) \ge 0$ and $K_h(z) \ge 0$ across all vertical levels.
2. **Energy Monotonicity:** Dissipation terms must satisfy non-negative column integration.
3. **Neutral Limit Recovery:** Closures must converge to standard MOST behavior as $Ri \to 0$.

---

### 3.3 $C^\infty$ Regularity & Smooth Operators

Numerical continuation engines (such as `BifurcationKit.jl`) require continuous Jacobians ($\mathbf{J} \in C^1$). ASM.jl replaces piecewise conditionals with algebraic $C^\infty$ smooth operators:

$$\operatorname{smooth\_max}(x, y; \epsilon) = \frac{x + y + \sqrt{(x - y)^2 + \epsilon^2}}{2}$$

Hyperbolic tangent weight functions blend stable and unstable similarity regimes smoothly:

$$w(\zeta) = \frac{1}{2} \left[ 1 + \tanh\left(\frac{\zeta}{\epsilon}\right) \right]$$

---

## 4. Discretization & Numerical Backends

ASM.jl implements two spatial discretization strategies:

```
                          [ Spatial Discretization ]
                                      │
         ┌────────────────────────────┴────────────────────────────┐
         ▼                                                         ▼
[ Stretched Finite Differences ]                         [ Spectral Gegenbauer-Galerkin ]
  • Physical grid with α stretching                        • Modal polynomial expansion C_n^(λ)
  • Direct Neumann boundary conditions                     • Precomputed 3-tensor C_ijk^(λ)
  • Stiff ODE scaling: λ_max ∝ Δz⁻²                        • Zero allocation in RHS calls (0 B)

```

### Spectral Galerkin Tensor Contractions

In `SpectralBLGalerkin`, non-linear variable-coefficient products are evaluated in modal space via precomputed 3-tensor contraction matrices:

$$C_{ijk}^{(\lambda)} = \int_{-1}^{1} C_i^{(\lambda)}(x) \, C_j^{(\lambda)}(x) \, C_k^{(\lambda)}(x) \, (1 - x^2)^{\lambda - 1/2} \, dx$$

This formulation eliminates dynamic memory allocations (`0 Bytes`) inside stiff time-stepping ODE loops (`Rodas5P`, `RadauIIA5`).

---

## 5. Model Selection, Calibration & Uncertainty Quantification

### 5.1 Information Criteria & Pareto Fronts

Candidate models identified across sparsity thresholds $\lambda$ are evaluated using information-theoretic criteria:

$$\text{AIC} = 2k + N \ln\left(\frac{\text{RSS}}{N}\right), \qquad \text{BIC} = k \ln(N) + N \ln\left(\frac{\text{RSS}}{N}\right)$$

Where:

* $k = \Vert{}\boldsymbol{\xi}\Vert{}_0$ is the active term count.
* $N$ is the number of weak-form projection rows.
* $\text{RSS} = \Vert{}\mathbf{b} - \mathbf{G}\boldsymbol{\xi}\Vert{}_2^2$ is the residual sum of squares.

`compute_pareto_front` extracts non-dominated models along the complexity ($k$) vs. error ($\text{RSS}$ or $R^2$) trade-off curve.

---

### 5.2 Hierarchical Bayesian Calibration

`HierarchicalTuring.jl` pools statistical strength across heterogeneous field campaigns:

```
                       [ Global Prior: θ_global ~ N(μ_0, Σ_0) ]
                                          │
         ┌────────────────────────────────┼────────────────────────────────┐
         ▼                                ▼                                ▼
  [ CASES-99 Site ]                  [ FLOSS Site ]                  [ SHEBA Site ]
θ_1 ~ N(θ_global, Σ_s)             θ_2 ~ N(θ_global, Σ_s)             θ_3 ~ N(θ_global, Σ_s)

```

Site parameters adapt to local surface roughness $z_0$ or soil moisture, while being constrained by the global population mean.

---

### 5.3 Profile Uncertainty Ribbons

`evaluate_profile_uncertainty` forwards parameter posterior chains $p(\boldsymbol{\theta} \mid \mathbf{y})$ into vertical eddy diffusivity distributions $K_m(z)$, returning $95\%$ credibility bands:

* **Median ($50.0\%$):** $Q_{0.50}\left(K_m(z \mid \boldsymbol{\theta}_i)\right)$
* **Lower Bound ($2.5\%$):** $Q_{0.025}\left(K_m(z \mid \boldsymbol{\theta}_i)\right)$
* **Upper Bound ($97.5\%$):** $Q_{0.975}\left(K_m(z \mid \boldsymbol{\theta}_i)\right)$

---

## 6. Observational Benchmark Diagnostics

ASM.jl has been validated against four field campaign datasets:

| Campaign | Observations ($N$) | Mean Stability ($\bar{\zeta}$) | Transversality ($\bar{\mathcal{T}}$) | Observed Boundary Layer Dynamics |
| --- | --- | --- | --- | --- |
| **FLOSS** | $70{,}796$ | $0.3842$ | $0.0068$ | Extreme manifold contraction over ice/snow; rapid mode damping. |
| **CASES-99** | $6{,}538$ | $0.2114$ | $0.1151$ | Strong non-equilibrium transients during sunset boundary layer collapse. |
| **BLLAST** | $5{,}600$ | $0.1985$ | $0.1053$ | Active fast-slow energy exchange across the evening transition. |
| **SHEBA** | $2{,}273$ | $0.3450$ | $0.0142$ | Strong manifold alignment in stable Arctic boundary layers. |

---

## 7. Developer & Workflow Guide

### 7.1 Key Makefile Targets

The `Makefile` automates execution, testing, and report generation:

```bash
# Display help and target menu
make help

# Resolve and instantiate dependencies
make instantiate

# Run standard test suite
make test

# Execute full pipeline with solver smoke tests
make test-smoke

# Run PDE closure benchmark & generate Pareto/LaTeX artifacts
make pde-benchmark

# Run multi-site campaign calibration & UQ export
make campaign-export

# Validate generated campaign tables and artifacts
make campaign-validate

# Clean generated reports
make clean

```

---

### 7.2 Automated Document Documentation Pipeline

To flatten `.md` files from `src/` into `docs/src/` with path attribution banners and collision-free names, run:

```bash
./scripts/copy_docs.sh

```

---

### 7.3 Executable Script Pipeline

To generate complete manuscript artifacts, run the three workflow scripts sequentially:

```bash
# 1. Benchmark discovery candidates, evaluate AIC/BIC & dump LaTeX term tables
julia --project=. scripts/run_pde_closure_benchmark.jl

# 2. Run MCMC/VI parameter calibration & output 95% K_m(z) credibility ribbons
julia --project=. scripts/run_campaign_exports.jl

# 3. Validate export schemas, files, and physical bounds
julia --project=. scripts/validate_campaign_exports.jl

```

All compiled equations (`.tex`), summary matrices (`.json`), and figure assets are saved directly into `reports/generated/`.