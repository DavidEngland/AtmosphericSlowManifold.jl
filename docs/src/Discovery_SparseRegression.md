<!-- Auto-generated from package source -->
> **Source:** `src/Discovery/SparseRegression.md`

`src/Discovery/SparseRegression.jl` implements sparse regression algorithms for model discovery in `AtmosphericSlowManifold.jl`, supporting both unconstrained iterative thresholding (STRidge) and physically constrained $L_1$-penalized Quadratic Programming (ConstrainedQP) solved via `JuMP.jl` and `HiGHS`.

---

### Sparse Optimization Algorithms

| Optimizer Struct | Optimization Paradigm | Penalty / Regularization | Physical Constraints Support |
| --- | --- | --- | --- |
| **`STRidge`** | Sequential Thresholded Ridge Regression | Iterative $L_2$ Ridge + hard thresholding ($\vert{}\xi_i\vert{} \ge \text{threshold}$) | No (Unconstrained) |
| **`ConstrainedQP`** | Convex Quadratic Program (JuMP + HiGHS) | $L_1$ LASSO penalty ($\lambda \Vert{}\boldsymbol{\xi}\Vert{}_1$) | Yes ($A_{\text{ineq}} \boldsymbol{\xi} \ge \mathbf{b}_{\text{ineq}}$) |

---

### Mathematical Formulations

#### 1. STRidge (Sequential Thresholded Ridge)

At iteration $k$, given current active support mask $S_k \subseteq \{1, \dots, n\}$, `STRidge` computes the ridge solution over active features:

$$\boldsymbol{\xi}_{S_k}^{(k)} = \left( \mathbf{G}_{:, S_k}^T \mathbf{G}_{:, S_k} + \lambda \mathbf{I} \right)^{-1} \mathbf{G}_{:, S_k}^T \mathbf{b}$$

Features failing the threshold test ($\vert{}\xi_i^{(k)}\vert{} < \text{threshold}$) are pruned from support set $S_{k+1}$. Iteration continues until support set stability ($S_{k+1} = S_k$) or $k = \text{max\_iter}$.

#### 2. ConstrainedQP (JuMP + HiGHS)

Solves the constrained LASSO optimization problem using auxiliary absolute value variables $t_i \ge \vert{}\xi_i\vert{}$:

$$\min_{\boldsymbol{\xi}, \mathbf{t}} \, \frac{1}{2} \Vert{}\mathbf{G}\boldsymbol{\xi} - \mathbf{b}\Vert{}_2^2 + \lambda \sum_{i=1}^{n_f} t_i$$

$$\text{subject to} \quad -\mathbf{t} \le \boldsymbol{\xi} \le \mathbf{t}, \qquad \mathbf{A}_{\text{ineq}} \boldsymbol{\xi} \ge \mathbf{b}_{\text{ineq}}$$

---

### API Reference & Method Dispatches

```julia
# 1. Unconstrained Sequential Thresholding
opt_stridge = STRidge(lambda = 1e-3, threshold = 1e-2, max_iter = 50)
xi_stridge = solve_sparse_regression(G, b, opt_stridge)

# 2. Physically Constrained Sparse QP
opt_qp = ConstrainedQP(lambda = 1e-3, solver = HiGHS.Optimizer)
xi_qp = solve_sparse_regression(G, b, opt_qp; A_ineq = A_mat, b_ineq = b_vec)

```