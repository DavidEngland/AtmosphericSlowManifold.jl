<!-- Auto-generated from package source -->
> **Source:** `src/Calibration/Backends/BayesianMCMC.md`

`src/Calibration/Backends/BayesianMCMC.jl` implements Bayesian parameter estimation and uncertainty quantification for discovered differential models, offering full MCMC posterior sampling via `Turing.jl` or a fast analytical Gaussian conjugate approximation fallback.

---

### Dual-Engine Inference Architecture

```
                       dispatch_calibrate(model, obs, alg::BayesianMCMC)
                                       │
                         Is Turing.jl available in env?
                                ┌──────┴──────┐
                             Yes│             │No / Fallback
                                ▼             ▼
                   _dispatch_calibrate_turing  _dispatch_calibrate_gaussian_approx
                   ┌────────────────────────┐  ┌─────────────────────────────────┐
                   │ • Dynamic Turing Model │  │ • Ridge Conjugate Posterior     │
                   │ • NUTS / HMC / MH      │  │ • Cholesky Sampling (F.L * z)   │
                   │ • Parameter Chain Draws│  │ • Exact Gaussian Covariance     │
                   └───────────┬────────────┘  └────────────────┬────────────────┘
                               └──────────────┬─────────────────┘
                                              ▼
                             CalibrationResult{BayesianMCMC}

```

---

### Mathematical Formulations (Gaussian Conjugate Fallback)

When `Turing.jl` is not installed or MCMC sampling fails, the engine falls back to analytical Gaussian linear-Bayesian regression over design matrix $\mathbf{X} \in \mathbb{R}^{N \times P}$ and observation target $\mathbf{y} \in \mathbb{R}^N$:

#### 1. Prior & Likelihood

* **Prior Distribution:** $\boldsymbol{\theta} \sim \mathcal{N}(\mathbf{0}, \lambda_{\text{prior}}^{-1} \mathbf{I}_p)$, with $\lambda_{\text{prior}} = 1.0$.
* **Likelihood:** $\mathbf{y} \mid \boldsymbol{\theta}, \sigma^2 \sim \mathcal{N}(\mathbf{X}\boldsymbol{\theta}, \sigma^2 \mathbf{I}_N)$.

#### 2. Posterior Distribution & Sampling

* **Posterior Precision Matrix:**

$$\mathbf{A} = \mathbf{X}^T \mathbf{X} + \lambda_{\text{prior}} \mathbf{I}_p$$


* **Posterior Mean ($\boldsymbol{\mu}_\theta$):**

$$\boldsymbol{\mu}_\theta = \mathbf{A}^{-1} \mathbf{X}^T \mathbf{y}$$


* **Residual Variance ($\sigma^2$):**

$$\sigma^2 = \max\left( \frac{\Vert{}\mathbf{y} - \mathbf{X}\boldsymbol{\mu}_\theta\Vert{}_2^2}{N}, 10^{-12} \right)$$


* **Posterior Covariance ($\boldsymbol{\Sigma}_\theta$):**

$$\boldsymbol{\Sigma}_\theta = \sigma^2 \mathbf{A}^{-1}$$


* **Posterior Parameter Draws:**
Given lower Cholesky factor $\mathbf{L}\mathbf{L}^T = \boldsymbol{\Sigma}_\theta + 10^{-10}\mathbf{I}$ and i.i.d. standard normal variates $\mathbf{z} \sim \mathcal{N}(\mathbf{0}, \mathbf{I})$:

$$\boldsymbol{\theta}^{(k)} = \boldsymbol{\mu}_\theta + \mathbf{L} \mathbf{z}^{(k)}$$



---

### Function Summary

| Function | Signature / Inputs | Operational Role |
| --- | --- | --- |
| **`_turing_available`** | `() -> Bool` | Checks whether `Turing.jl` is available in the current Julia load path. |
| **`_bayes_draws`** | `(mu, cov, ndraws)` | Draws $K$ posterior parameter samples via Cholesky decomposition of symmetric covariance matrix $\boldsymbol{\Sigma}_\theta$. |
| **`_rowwise_std`** | `(draws) -> Vector{Float64}` | Computes empirical sample standard deviations across generated parameter MCMC/analytical draws. |
| **`_extract_theta_draws`** | `(chain, p)` | Extracts per-parameter draw vectors (`c_1, c_2, ...`) from a `Turing.Chains` object. |
| **`_dispatch_calibrate_turing`** | `(model, obs, alg)` | Evaluates a runtime Turing probabilistic model (`bayes_wsindy_kernel`) and samples using `NUTS`, `HMC`, or `MH`. |
| **`_dispatch_calibrate_gaussian_approx`** | `(model, obs, alg)` | Computes analytical conjugate posterior mean, variance, and Cholesky draws when `Turing.jl` is absent. |
| **`dispatch_calibrate`** | `(model, obs, alg::BayesianMCMC)` | Main dispatch entry point routing execution to Turing or Gaussian approximation. |

---

### Result Payload & Parameter Dictionary

The function returns a `CalibrationResult{BayesianMCMC}` containing:

* **`parameters` Dictionary:**
* `:posterior_mean` — $P$-element vector of estimated mean coefficients.
* `:posterior_std` — $P$-element vector of standard errors / uncertainties.
* `:c_1, :c_2, ...` — Individual sample chains across all iterations ($K = \text{samples} \times \text{chains}$).


* **`diagnostics` Dictionary:**
* Engine provenance (`:turing` vs. `:gaussian_approx`).
* Sampler specifications (`:NUTS`, `:HMC`, or `:MH`, `target_accept`, `samples`, `chains`).
* Residual error statistics and sample metadata.# s