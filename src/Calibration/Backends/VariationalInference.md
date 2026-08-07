`src/Calibration/Backends/VariationalInference.jl` implements fast variational Bayesian parameter estimation for discovered models. It extends `dispatch_calibrate` for the `VariationalInference` algorithm type, optimizing an Evidence Lower Bound (ELBO) over a mean-field Gaussian variational family $q(\boldsymbol{\theta})$ using ADAM or SGD gradient descent.

---

### Mathematical Foundation & Variational Mechanics

#### 1. Mean-Field Variational Family & Prior Form

The posterior distribution $p(\boldsymbol{\theta} \mid \mathbf{y})$ over coefficient vector $\boldsymbol{\theta} \in \mathbb{R}^P$ is approximated by a fully factorized (mean-field) multivariate Gaussian:

$$q(\boldsymbol{\theta}) = \mathcal{N}(\boldsymbol{\mu}, \operatorname{diag}(\boldsymbol{\sigma}^2)) = \prod_{j=1}^P \mathcal{N}(\theta_j; \mu_j, \sigma_j^2)$$

* **Prior Distribution:** Isotropic Gaussian prior $p(\boldsymbol{\theta}) = \mathcal{N}(\mathbf{0}, \tau^2 \mathbf{I})$ with default variance scale $\tau^2 = 1.0$.
* **Likelihood Noise:** Normal noise model $\mathbf{y} \mid \mathbf{X}, \boldsymbol{\theta} \sim \mathcal{N}(\mathbf{X}\boldsymbol{\theta}, \sigma_y^2 \mathbf{I})$, where $\sigma_y^2$ is initialized from the sample variance of target observations $\mathbf{y}$.

#### 2. Analytical Evidence Lower Bound (ELBO)

The ELBO $\mathcal{L}(\boldsymbol{\mu}, \boldsymbol{\sigma})$ decouples into expected log-likelihood and analytical Kullback–Leibler (KL) divergence terms:

$$\mathcal{L}(\boldsymbol{\mu}, \boldsymbol{\sigma}) = \mathbb{E}_{q(\boldsymbol{\theta})}[\log p(\mathbf{y} \mid \mathbf{X}, \boldsymbol{\theta})] - D_{\text{KL}}(q(\boldsymbol{\theta}) \parallel p(\boldsymbol{\theta}))$$

1. **Expected Squared Error:**

$$\mathbb{E}_q \left[ \Vert{}\mathbf{y} - \mathbf{X}\boldsymbol{\theta}\Vert{}_2^2 \right] = \Vert{}\mathbf{y} - \mathbf{X}\boldsymbol{\mu}\Vert{}_2^2 + \sum_{j=1}^P g_j \sigma_j^2, \qquad g_j = \sum_{i=1}^N X_{i, j}^2$$


2. **Kullback–Leibler Divergence:**

$$D_{\text{KL}}(q \parallel p) = \frac{1}{2} \sum_{j=1}^P \left[ \frac{\sigma_j^2 + \mu_j^2}{\tau^2} - 1 - \log \left( \frac{\sigma_j^2}{\tau^2} \right) \right]$$



Putting it together:

$$\mathcal{L}(\boldsymbol{\mu}, \boldsymbol{\sigma}) = -\frac{1}{2 \sigma_y^2} \left( \Vert{}\mathbf{y} - \mathbf{X}\boldsymbol{\mu}\Vert{}_2^2 + \sum_{j=1}^P g_j \sigma_j^2 \right) - \frac{1}{2} \sum_{j=1}^P \left[ \frac{\sigma_j^2 + \mu_j^2}{\tau^2} - 1 - \log \left( \frac{\sigma_j^2}{\tau^2} \right) \right]$$

#### 3. Exact Gradient Formulation

To minimize negative ELBO, exact analytical gradients are derived with respect to mean parameters $\boldsymbol{\mu}$ and log-scale parameters $s_j = \log \sigma_j$:

$$\mathbf{g}_{\boldsymbol{\mu}} = -\frac{\mathbf{X}^T (\mathbf{y} - \mathbf{X}\boldsymbol{\mu})}{\sigma_y^2} + \frac{\boldsymbol{\mu}}{\tau^2}$$

$$\mathbf{g}_{\mathbf{s}} = \frac{\mathbf{g} \odot \boldsymbol{\sigma}^2}{\sigma_y^2} + \frac{\boldsymbol{\sigma}^2}{\tau^2} - \mathbf{1}$$

---

### Algorithm Architecture & Calibration Flow

```
                                  dispatch_calibrate(model, obs, alg::VariationalInference)
                                                              │
                                            Construct Design Matrix X and Target y
                                                              │
                                            Initialize μ via Ridge OLS & log(σ) = -2.30
                                                              │
                                            ┌─────────────────┴─────────────────┐
                                            │ Optimization Loop (ADAM / SGD)    │
                                            │                                   │
                                            │ 1. Evaluate Analytical ELBO       │
                                            │ 2. Compute Gradients (g_μ, g_s)   │
                                            │ 3. Check ||g||_2 < 1e-6           │
                                            │ 4. Update μ and log(σ)            │
                                            │ 5. Clamp log(σ) ∈ [-13.82, 2.30]  │
                                            └─────────────────┬─────────────────┘
                                                              │
                                            Generate N_draws Monte Carlo Samples
                                                              │
                                           CalibrationResult{VariationalInference}

```

---

### Calibration Result Output Specs

| Component | Key / Type | Mathematical / Operational Content |
| --- | --- | --- |
| **`model`** | `DiscoveredModel{Float64}` | Updated model instance containing variational mean parameters $\boldsymbol{\mu}$ as active coefficients. |
| **`parameters[:variational_mean]`** | `Vector{Float64}` | Estimated posterior mean vector $\boldsymbol{\mu} \in \mathbb{R}^P$. |
| **`parameters[:variational_scale]`** | `Vector{Float64}` | Estimated posterior standard deviation vector $\boldsymbol{\sigma} \in \mathbb{R}^P$. |
| **`parameters[:c_1, :c_2, ...]`** | `Vector{Float64}` | $N_{\text{samples}}$ Monte Carlo draws generated from $q(\theta_i) = \mathcal{N}(\mu_i, \sigma_i^2)$ for term $i$. |
| **`diagnostics[:elbo_trace]`** | `Vector{Float64}` | Iteration-by-iteration history of the ELBO objective to verify optimization convergence. |
| **`diagnostics[:converged]`** | `Bool` | Convergence flag indicating whether gradient norm fell below threshold $10^{-6}$. |
| **`diagnostics[:noise_var]`** | `Float64` | Target observation noise variance $\sigma_y^2$ used during optimization. |