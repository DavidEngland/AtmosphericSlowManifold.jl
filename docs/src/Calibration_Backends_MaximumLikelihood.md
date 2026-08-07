<!-- Auto-generated from package source -->
> **Source:** `src/Calibration/Backends/MaximumLikelihood.md`

`src/Calibration/Backends/MaximumLikelihood.jl` implements deterministic point estimation for discovered PDE model coefficients via Ridge-regularized Maximum Likelihood / Ordinary Least Squares regression. It extends `dispatch_calibrate` for the `MaximumLikelihood` algorithm type.

---

### Mathematical Foundation & Estimation Pipeline

#### 1. Design Matrix Assembly ($\mathbf{X}$)

Given an `ObservationTable` and a `DiscoveredModel` containing terms $T_j = \prod_k f_{j,k}^{p_{j,k}}$, the backend constructs a design matrix $\mathbf{X} \in \mathbb{R}^{N \times P}$ where column $j$ corresponds to evaluated basis product values across all $N$ observation samples:

$$\mathbf{X}_{i, j} = \prod_{b \in T_j.\text{basis}} f_b(z_i, t_i)^{p_b}$$

#### 2. Regularized Maximum Likelihood Solves ($\hat{\boldsymbol{\theta}}_{\text{MLE}}$)

To prevent ill-conditioned system inversions when terms in candidate libraries exhibit collinearity, parameter vectors are estimated using $L_2$-regularized Ridge estimation:

$$\hat{\boldsymbol{\theta}}_{\text{MLE}} = \left( \mathbf{X}^T \mathbf{X} + \lambda \mathbf{I} \right)^{-1} \mathbf{X}^T \mathbf{y}$$

where regularization scalar $\lambda = \max(\text{alg.tol}, 10^{-12})$.

#### 3. Diagnostic & Goodness-of-Fit Computation

After estimating parameters, the module computes residual sum of squares ($\text{RSS}$) and coefficient of determination ($R^2$):

$$\text{RSS} = \Vert{}\mathbf{y} - \mathbf{X}\hat{\boldsymbol{\theta}}_{\text{MLE}}\Vert{}_2^2, \qquad R^2 = 1 - \frac{\text{RSS}}{\sum_{i=1}^N (y_i - \bar{y})^2}$$

---

### Function Summary

| Function | Signature / Inputs | Operational Role |
| --- | --- | --- |
| **`_cal_feature_column`** | `(obs, feature)` | Extracts observation column data or evaluates vertical spatial derivatives ($\partial_z^k f$) via finite differences. |
| **`_cal_term_column`** | `(obs, term)` | Evaluates composite product term $T_j = \prod f_b^{p_b}$ vector across observation nodes. |
| **`_cal_design_matrix`** | `(obs, model)` | Assembles full feature matrix $\mathbf{X} \in \mathbb{R}^{N \times P}$ for all model terms. |
| **`_cal_fit_ridge`** | `(X, y, λ)` | Solves regularized normal equations $(\mathbf{X}^T \mathbf{X} + \lambda \mathbf{I}) \hat{\boldsymbol{\theta}} = \mathbf{X}^T \mathbf{y}$. |
| **`_cal_update_model`** | `(model, θ, residual_norm)` | Constructs updated `DiscoveredModel` holding calibrated coefficients $\boldsymbol{\theta}$. |
| **`dispatch_calibrate`** | `(model, obs, alg::MaximumLikelihood)` | Executes point estimation pipeline and returns `CalibrationResult{MaximumLikelihood}` with $R^2$ diagnostics. |