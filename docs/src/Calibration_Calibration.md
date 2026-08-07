<!-- Auto-generated from package source -->
> **Source:** `src/Calibration/Calibration.md`

`src/Calibration/Calibration.jl` defines the parameter estimation and uncertainty quantification submodule for `AtmosphericSlowManifold.jl`, re-exporting symbolic features and discovered model representations to calibrate discovered PDE coefficients against observational data.

---

### Submodule Architecture

| Inclusion File | Subsystem Functionality | Estimation & Mathematical Role |
| --- | --- | --- |
| **`Interface.jl`** | Common API & Cost Functions | Defines standard calibration problem structures, parameter bounds, residual evaluation protocols, and likelihood functions $p(\mathbf{y} \mid \boldsymbol{\theta})$. |
| **`Backends/MaximumLikelihood.jl`** | Point Estimation | Fits deterministic optimal parameter values $\hat{\boldsymbol{\theta}}_{\text{MLE}} = \arg\max_{\boldsymbol{\theta}} \log p(\mathbf{y} \mid \boldsymbol{\theta})$ via gradient-based or derivative-free optimization. |
| **`Backends/BayesianMCMC.jl`** | Posterior Sampling | Sample-based uncertainty quantification computing full parameter posterior distributions $p(\boldsymbol{\theta} \mid \mathbf{y}) \propto p(\mathbf{y} \mid \boldsymbol{\theta}) p(\boldsymbol{\theta})$ using MCMC. |
| **`Backends/VariationalInference.jl`** | Approximate Bayesian Inference | Fast variational approximation fitting tractable distributions $q(\boldsymbol{\theta}) \approx p(\boldsymbol{\theta} \mid \mathbf{y})$ by maximizing the Evidence Lower Bound (ELBO). |

---

### Calibration Pipeline Workflow

```
[ Discovery Engine ]   ──> DiscoveredModel ──┐
                                             ├──> [ Calibration.jl ] ──> Interface.jl
[ Observation Data ]   ──> ObservationTable ──┘                             │
                                                                            ├──> MaximumLikelihood.jl   ──> θ_MLE
                                                                            ├──> BayesianMCMC.jl        ──> Posterior Samples
                                                                            └──> VariationalInference.jl──> q(θ) Density

```

Would you like to build out `src/Calibration/Interface.jl` or implement one of the estimation backends next?