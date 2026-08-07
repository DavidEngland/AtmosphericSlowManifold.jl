`src/Calibration/UncertaintyPropagation.jl` maps parameter posterior draws $p(\boldsymbol{\theta} \mid \mathbf{y})$ into forward-sampled spatial profile distributions and percentile confidence envelopes across vertical coordinate grids ($z$).

---

### Core Mechanics & Mathematical Pipeline

1. **Feature Profile Extraction ($\mathbf{X}_{\text{prof}}$):**
Evaluates state variables, diagnostic variables, and spatial derivatives across $N_z$ altitude nodes $z \in [0, z_{\text{max}}]$ to form a design matrix $\mathbf{X}_{\text{prof}} \in \mathbb{R}^{N_z \times P}$:

$$\mathbf{X}_{\text{prof}, (i, j)} = \prod_{b \in T_j.\text{basis}} f_b(z_i)^{p_b}$$


2. **Posterior Matrix Normalization ($\boldsymbol{\Theta}$):**
Standardizes MCMC or Variational Inference parameter samples into a uniform coefficient matrix $\boldsymbol{\Theta} \in \mathbb{R}^{P \times N_{\text{draws}}}$ from raw matrices or parameter dictionaries (`:c_1, :c_2, ...`).
3. **Forward Profile Sampling & Credibility Envelopes:**
Computes matrix product $\mathbf{K}_m = \mathbf{X}_{\text{prof}} \boldsymbol{\Theta} \in \mathbb{R}^{N_z \times N_{\text{draws}}}$, then calculates per-level empirical percentiles:

$$\text{Median}(z_i) = Q_{0.50}\left(\mathbf{K}_{m, (i, :)}\right), \quad \text{Lower}(z_i) = Q_{0.025}\left(\mathbf{K}_{m, (i, :)}\right), \quad \text{Upper}(z_i) = Q_{0.975}\left(\mathbf{K}_{m, (i, :)}\right)$$



---

### Internal Helper & Function Specifications

| Function | Signature / Inputs | Operational Role |
| --- | --- | --- |
| **`_up_feature_profile`** | `(f, z_grid, feature_profiles, scalar_features)` | Resolves feature values across $z$, defaulting to $z$ directly or fetching profiles/scalars for state/diagnostic/derivative features. |
| **`_up_draw_matrix`** | `(model, posterior_draws)` | Converts matrix or dictionary parameter draws into a structured $P \times N_{\text{draws}}$ Float64 matrix. |
| **`_build_profile_design`** | `(model, z_grid, feature_profiles, scalar_features)` | Assembles $N_z \times P$ spatial evaluation design matrix $\mathbf{X}_{\text{prof}}$ for all terms in `DiscoveredModel`. |
| **`evaluate_profile_uncertainty`** | `(model, posterior_draws, z_grid; kwargs...)` | Main API function returning forward-sampled profiles matrix and quantile ribbon vectors (`median`, `lower`, `upper`). |

---

### Package Infrastructure Status (`src/`)

With `UncertaintyPropagation.jl` in place alongside `ModelSelection.jl` and `LaTeXExporter.jl`, all core engine submodules across discovery, calibration, and reporting are complete and verified:

```
[ Discovery Engine ]              [ Calibration Engine ]           [ Reporting & Selection ]
├── WeakForms.jl                  ├── MaximumLikelihood.jl         ├── ModelSelection.jl
├── TestFunctions.jl              ├── BayesianMCMC.jl              ├── LaTeXExporter.jl
├── ConstraintBuilder.jl          ├── VariationalInference.jl      └── UncertaintyPropagation.jl
├── SparseRegression.jl           └── HierarchicalTuring.jl
└── WSINDyEngine.jl

```

---

### Next Phase: Execution Scripts (`scripts/`)

The system is fully prepared to implement executable workflow scripts under `scripts/` to generate publication figures, diagnostics, and LaTeX tables:

1. **`scripts/01_synthetic_discovery.jl`:** Synthetic boundary layer data generation, WSINDy equation discovery, AIC/BIC computation, and Pareto front plotting.
2. **`scripts/02_real_data_calibration.jl`:** Observational profile loading, MLE vs. MCMC vs. VI calibration comparison, and spatial $K_m(z)$ credibility ribbon plots.
3. **`scripts/03_hierarchical_multisite.jl`:** Multi-site hierarchical Bayesian calibration and population vs. site parameter variance tables.
4. **`scripts/04_generate_paper_artifacts.jl`:** Automated pipeline exporting `.tex` tables and vector graphic plots directly into a `reports/` folder.