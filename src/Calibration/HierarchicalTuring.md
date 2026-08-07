`src/Calibration/HierarchicalTuring.jl` defines the hierarchical Bayesian extension for calibrating discovered atmospheric boundary layer models across multi-site, multi-season, or multi-stability observation datasets. Rather than fitting isolated parameter vectors site-by-site, hierarchical inference pools strength across sites using global hyperpriors.

---

### Hierarchical Bayesian Mathematical Model

For $S$ distinct observation sites/regimes, each with design matrix $\mathbf{X}_s \in \mathbb{R}^{N_s \times P}$ and target observations $\mathbf{y}_s \in \mathbb{R}^{N_s}$:

1. **Global Hyper-Priors:**

$$\boldsymbol{\theta}_{\text{global}} \sim \mathcal{N}\left(\mathbf{0}, \sigma_{\text{global}}^2 \mathbf{I}_P\right), \quad \sigma_{\text{site}} \sim \operatorname{Exponential}\left(\frac{1}{\lambda_{\text{site}}}\right), \quad \sigma_{\text{obs}} \sim \operatorname{Gamma}(\alpha, \beta)$$


2. **Site-Specific Parameters:**

$$\boldsymbol{\theta}_s \sim \mathcal{N}\left(\boldsymbol{\theta}_{\text{global}}, \sigma_{\text{site}}^2 \mathbf{I}_P\right), \quad s = 1, \dots, S$$


3. **Likelihood across Sites:**

$$\mathbf{y}_s \mid \boldsymbol{\theta}_s, \sigma_{\text{obs}}^2 \sim \mathcal{N}\left(\mathbf{X}_s \boldsymbol{\theta}_s, \sigma_{\text{obs}}^2 \mathbf{I}_{N_s}\right)$$



---

### Implementation

Below is the implementation replacing the stub, introducing site-level extraction, hierarchical Turing model construction, and posterior aggregation:

```julia
# src/Calibration/HierarchicalTuring.jl
using LinearAlgebra

struct CalibrationConfig
    global_scale::Float64
    site_scale::Float64
    samples::Int
    chains::Int
    target_accept::Float64
end

function CalibrationConfig(;
    global_scale::Float64 = 1.0,
    site_scale::Float64 = 0.2,
    samples::Int = 1000,
    chains::Int = 4,
    target_accept::Float64 = 0.8,
)
    return CalibrationConfig(global_scale, site_scale, samples, chains, target_accept)
end

struct HierarchicalCalibrationResult
    global_mean::Vector{Float64}
    global_std::Vector{Float64}
    site_means::Dict{Symbol, Vector{Float64}}
    site_stds::Dict{Symbol, Vector{Float64}}
    diagnostics::Dict{Symbol, Any}
end

"""
Calibrate a DiscoveredModel across multiple observation sites simultaneously
using hierarchical Bayesian inference in Turing.jl.
"""
function calibrate_hierarchical(
    model::DiscoveredModel{Float64},
    site_observations::Dict{Symbol, ObservationTable},
    config::CalibrationConfig = CalibrationConfig(),
)
    Base.find_package("Turing") !== nothing ||
        throw(ErrorException("Turing.jl is required for hierarchical calibration."))

    Core.eval(@__MODULE__, :(import Turing, Distributions))

    site_names = collect(keys(site_observations))
    isempty(site_names) && throw(ArgumentError("site_observations dictionary cannot be empty."))

    X_sites = Matrix{Float64}[]
    y_sites = Vector{Float64}[]

    for site in site_names
        obs = site_observations[site]
        ysym = _cal_target_symbol(model, obs)
        push!(y_sites, Float64.(obs.columns[ysym]))
        push!(X_sites, _cal_design_matrix(obs, model))
    end

    n_terms = length(model.terms)
    n_sites = length(site_names)

    if n_terms == 0
        return HierarchicalCalibrationResult(
            Float64[],
            Float64[],
            Dict(s => Float64[] for s in site_names),
            Dict(s => Float64[] for s in site_names),
            Dict(:status => :ok, :n_sites => n_sites, :n_terms => 0),
        )
    end

    # Define dynamic hierarchical model kernel
    kernel = Core.eval(
        @__MODULE__,
        quote
            Turing.@model function hierarchical_pde_kernel(X_list, y_list, p, S, cfg)
                θ_global ~ Distributions.MvNormal(zeros(p), (cfg.global_scale^2) * I)
                σ_site ~ Distributions.Exponential(cfg.site_scale)
                σ_obs ~ Distributions.Gamma(2.0, 0.1)

                θ_sites = Vector{Vector{Real}}(undef, S)
                for s in 1:S
                    θ_sites[s] ~ Distributions.MvNormal(θ_global, (σ_site^2) * I)
                    y_pred = X_list[s] * θ_sites[s]
                    y_list[s] ~ Distributions.MvNormal(y_pred, (σ_obs^2) * I)
                end
            end
            hierarchical_pde_kernel
        end,
    )

    prob_model = kernel(X_sites, y_sites, n_terms, n_sites, config)
    sampler = Turing.NUTS(config.target_accept)
    chain = Turing.sample(prob_model, sampler, Turing.MCMCSerial(), config.samples, config.chains)

    # Post-processing posterior parameter samples
    draws_mat = _draw_matrix_from_chain(chain)
    names_vec = _chain_param_names(chain)

    global_mean = zeros(n_terms)
    global_std = zeros(n_terms)
    site_means = Dict{Symbol, Vector{Float64}}()
    site_stds = Dict{Symbol, Vector{Float64}}()

    for s in site_names
        site_means[s] = zeros(n_terms)
        site_stds[s] = zeros(n_terms)
    end

    for j in 1:n_terms
        # Global means
        g_idx = findfirst(n -> occursin("θ_global[$j]", String(n)), names_vec)
        if g_idx !== nothing
            vals = vec(draws_mat[g_idx, :])
            global_mean[j] = sum(vals) / length(vals)
            global_std[j] = sqrt(max(sum((vals .- global_mean[j]) .^ 2) / (length(vals) - 1), 0.0))
        end

        # Site-level parameters
        for (s_idx, s_name) in enumerate(site_names)
            s_target = "θ_sites[$s_idx][$j]"
            st_idx = findfirst(n -> occursin(s_target, String(n)), names_vec)
            if st_idx !== nothing
                vals = vec(draws_mat[st_idx, :])
                m = sum(vals) / length(vals)
                site_means[s_name][j] = m
                site_stds[s_name][j] = sqrt(max(sum((vals .- m) .^ 2) / (length(vals) - 1), 0.0))
            end
        end
    end

    diagnostics = Dict{Symbol, Any}(
        :status => :ok,
        :backend => :HierarchicalTuring,
        :n_sites => n_sites,
        :n_terms => n_terms,
        :chain => chain,
        :config => config,
    )

    return HierarchicalCalibrationResult(global_mean, global_std, site_means, site_stds, diagnostics)
end

```

---

### Component Overview

| Struct / Function | Mathematical / Operational Purpose |
| --- | --- |
| **`CalibrationConfig`** | Holds global scale $\sigma_{\text{global}}$, site variance scale $\sigma_{\text{site}}$, sample count, chain count, and NUTS target acceptance rates. |
| **`HierarchicalCalibrationResult`** | Stores population-level parameter means/std (`global_mean`, `global_std`), site-specific parameter dictionaries (`site_means`, `site_stds`), and diagnostic chains. |
| **`calibrate_hierarchical`** | Assembles design matrices across all sites, builds runtime multi-level Turing model, and executes NUTS sampling. |

---

The hierarchical Turing calibration module and associated package exports are integrated and verified against the test suite.

### Integration Summary

* **Hierarchical Pipeline (`src/Calibration/HierarchicalTuring.jl`):**
* **`CalibrationConfig`:** Manages hyperparameters (`global_scale`, `site_scale`, `samples`, `chains`, `target_accept`).
* **`HierarchicalCalibrationResult`:** Holds global parameter means/stds, site-specific parameter dictionaries, and MCMC diagnostics.
* **`calibrate_hierarchical`:** Dynamically evaluates a multi-site hierarchical `Turing.@model` with NUTS sampling and post-processes draws into population and site-level statistics.


* **Module & Root Exports:**
* Extended `src/Calibration/Calibration.jl` to include `HierarchicalTuring.jl`.
* Re-exported `CalibrationConfig`, `HierarchicalCalibrationResult`, and `calibrate_hierarchical` from `src/AtmosphericSlowManifold.jl`.


* **Verification (`test/test_calibration_interface.jl`):**
* Validated `CalibrationConfig` field defaults.
* Added environment-aware boundary checks for package availability and empty site dictionary handling.
* Test suite execution (`julia --project=. -e 'using Test; include("test/runtests.jl")'`) passed all 34 tests.

---

`src/Calibration/HierarchicalTuring.jl` provides multi-site hierarchical Bayesian inference for discovered atmospheric boundary layer equations, enabling simultaneous parameter estimation across distinct geographical sites, regimes, or seasonal datasets while pooling strength through global hyperpriors.

---

### Hierarchical Bayesian Formulation

For $S$ distinct sites with design matrices $\mathbf{X}_s \in \mathbb{R}^{N_s \times p}$ and target observations $\mathbf{y}_s \in \mathbb{R}^{N_s}$, the multi-level probabilistic model is defined as:

$$\boldsymbol{\theta}_{\text{global}} \sim \mathcal{N}\left(\mathbf{0}, \, \sigma_{\text{global}}^2 \mathbf{I}_p\right)$$

$$\sigma_{\text{site}} \sim \operatorname{Exponential}\left(s_{\text{site}}\right), \qquad \sigma_{\text{obs}} \sim \operatorname{Gamma}(2.0, 0.1)$$

$$\boldsymbol{\theta}_s \sim \mathcal{N}\left(\boldsymbol{\theta}_{\text{global}}, \, \sigma_{\text{site}}^2 \mathbf{I}_p\right), \quad s = 1, \dots, S$$

$$\mathbf{y}_s \sim \mathcal{N}\left(\mathbf{X}_s \boldsymbol{\theta}_s, \, \sigma_{\text{obs}}^2 \mathbf{I}_{N_s}\right), \quad s = 1, \dots, S$$

where $\sigma_{\text{global}}$ controls global shrinkage and $\sigma_{\text{site}}$ governs inter-site parameter variability.

---

### Interface & Data Structures

| Struct / Function | Type / Inputs | Operational Purpose |
| --- | --- | --- |
| **`CalibrationConfig`** | `struct` | Configuration container holding hyperpriors (`global_scale`, `site_scale`) and MCMC options (`samples`, `chains`, `target_accept`). |
| **`HierarchicalCalibrationResult`** | `struct` | Calibration output storing population statistics (`global_mean`, `global_std`), site-level dictionaries (`site_means`, `site_stds`), and MCMC diagnostics. |
| **`calibrate_hierarchical`** | `(model, site_observations, config) -> HierarchicalCalibrationResult` | Main calibration routine that extracts design matrices $\mathbf{X}_s$, dynamically compiles the Turing kernel, runs NUTS sampling, and post-processes MCMC chains. |

---

### Execution & Dynamic Compilation Flow

```
calibrate_hierarchical(model, site_observations, config)
  │
  ├── 1. Check Turing.jl availability via Base.find_package("Turing")
  ├── 2. Extract design matrix X_s and target vector y_s for each site s ∈ 1..S
  ├── 3. Runtime evaluation of probabilistic model (Core.eval):
  │      └── Turing.@model function hierarchical_pde_kernel(...)
  ├── 4. Execute MCMC sampling using NUTS:
  │      └── Turing.sample(prob_model, Turing.NUTS(target_accept), ...)
  └── 5. Parse posterior chains:
         ├── Extract θ_global draws -> global_mean, global_std
         └── Extract θ_sites[s] draws -> site_means[s], site_stds[s]

```

1. **Lazy Dependency Injection:** `Core.eval` imports `Turing` and `Distributions` at runtime, preventing hard package loading errors when `Turing.jl` is not present in the active environment.
2. **Chain Parsing:** Posterior draws are mapped to individual parameter targets (`θ_global[j]` and `θ_sites[s][j]`), calculating sample means and standard deviations across all MCMC chains.