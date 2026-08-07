<!-- Auto-generated from package source -->
> **Source:** `src/Discovery/ConstraintBuilder.md`

`src/Discovery/ConstraintBuilder.jl` constructs physical linear inequality constraints ($A_{\text{ineq}} \boldsymbol{\xi} \ge \mathbf{b}_{\text{ineq}}$) to ensure discovered SINDy models obey fundamental thermodynamic principles, non-negative dissipation, and monotonicity bounds.

---

### Physical Constraint Taxonomy

| Constraint Type | Struct Definition | Mathematical / Physical Condition | Linear Form ($A_{k,:} \boldsymbol{\xi} \ge b_k$) |
| --- | --- | --- | --- |
| **Positivity** | `PositivityConstraint(target)` | Requires target coefficient $\xi_i \ge 0$ (e.g., non-negative eddy diffusivities $K_m, K_h$). | $\xi_i \ge 0$ |
| **Monotonicity** | `MonotonicityConstraint(target, respect_to)` | Enforces sign agreement of target coefficient $\xi_i$ relative to gradient orientation $\text{sgn}\left(\sum_k f_{j,k}\right)$. | $\operatorname{sgn}\left(\bar{f}_j\right) \xi_i \ge 0$ |
| **Energy** | `EnergyConstraint(dissipation_term)` | Restricts turbulent energy dissipation terms $\xi_i \ge 0$ to prevent unphysical energy creation. | $\xi_i \ge 0$ |

---

### Constraint Matrix Assembly Flow

1. **Feature Identification (`_feature_index`):**
Searches the candidate feature library for exact type and structural equality (`typeof(f) == typeof(target) && f == target`), returning matrix column index $i$.
2. **System Matrix Construction (`assemble_constraint_matrix`):**
* Maps each constraint into an inequality coefficient row vector $\mathbf{a}_k \in \mathbb{R}^{1 \times N_{\text{features}}}$ and corresponding right-hand-side scalar $b_k = 0.0$.
* Concatenates row vectors into `PhysicalConstraintMatrix`:

$$A_{\text{ineq}} \boldsymbol{\xi} \ge \mathbf{b}_{\text{ineq}}$$




3. **Empty Fallback Handling:**
If no constraints are specified, returns an empty matrix `zeros(0, N_features)` and empty vector `Float64[]`, allowing downstream sparse regression solvers (e.g., SR3) to bypass inequality projections.