`src/Calibration/Interface.jl` defines the dispatch contracts, algorithm configuration types, and result containers for calibrating `DiscoveredModel` parameters against `ObservationTable` datasets.

---

### Calibration Algorithm Types & Data Structures

| Struct / Interface | Supertype / Fields | Operational Purpose |
| --- | --- | --- |
| **`AbstractCalibrationAlgorithm`** | Abstract Base | Parent type for point estimation, sample-based Bayesian inference, and variational approximation algorithms. |
| **`BayesianMCMC`** | `<: AbstractCalibrationAlgorithm` | MCMC configuration (`samples`, `chains`, `target_accept`, `sampler_type`, e.g., `:NUTS`). |
| **`MaximumLikelihood`** | `<: AbstractCalibrationAlgorithm` | Point estimation setup (`optimizer`, `max_iter`, `tol`). |
| **`VariationalInference`** | `<: AbstractCalibrationAlgorithm` | Variational Bayes setup (`samples`, `max_iter`, `optimizer`). |
| **`CalibrationResult{A}`** | Parametric Struct | Holds calibrated `model::DiscoveredModel{Float64}`, estimated parameter arrays (`parameters`), convergence/MCMC metrics (`diagnostics`), and algorithm instance (`algorithm::A`). |

---

### Dispatch Pipeline & Call Flow

```
calibrate(model, obs; algorithm = BayesianMCMC(...))
  │
  └──> dispatch_calibrate(model, obs, algorithm)
         ├──> Backends/MaximumLikelihood.jl (MLE / MAP)
         ├──> Backends/BayesianMCMC.jl (NUTS / HMC)
         └──> Backends/VariationalInference.jl (ADVI / ELBO)

```

1. **User Entry Point (`calibrate`):** Provides a unified user API defaulting to NUTS MCMC sampling (`samples = 1000`, `chains = 4`, `target_accept = 0.8`).
2. **Backend Extension Point (`dispatch_calibrate`):** Throws a default `error(...)` until concrete backend implementation files (`MaximumLikelihood.jl`, `BayesianMCMC.jl`, `VariationalInference.jl`) extend `dispatch_calibrate` for their respective algorithm types.

Would you like to implement one of the backend modules next, such as `src/Calibration/Backends/MaximumLikelihood.jl` or `src/Calibration/Backends/BayesianMCMC.jl`?