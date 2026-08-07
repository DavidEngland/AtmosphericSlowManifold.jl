`src/Discovery/WSINDyEngine.jl` orchestrates the end-to-end Weak-form Sparse Identification of Nonlinear Dynamics (WSINDy) pipeline, integrating weak-form integration, physical constraint matrices, sparse optimization, and symbolic model/closure construction into a unified execution layer.

---

### End-to-End Discovery Workflows

The module supports two distinct discovery workflows:

```
[ ObservationTable ] ──> [ FeatureLibrary / Candidates ]
                               │
               ┌───────────────┴───────────────┐
               ▼                               ▼
       discover(...)                   discover_closure(...)
  (Typed Feature-Library Path)    (Symbolic Vector{Num} Path)
               │                               │
       Weak Assembly                   Weak Assembly
               │                               │
       Constraint Matrix               Positivity / Custom Bounds
               │                               │
      STRidge / ConstrainedQP             JuMP / HiGHS (QP)
               │                               │
               ▼                               ▼
       DiscoveredModel                   WSINDyClosure

```

---

### Core Pipeline Interface

#### 1. Feature Library Discovery (`discover`)

Discovers general governing differential equations or diagnostic representations from structured `FeatureLibrary` objects:

```julia
model = discover(
    obs,
    library,
    constraints,
    test_family,
    optimizer;
    target_variable = :K_m,
    target = :u,
    threshold = 1e-8
)

```

1. **Weak-Form Assembly:** Integrates spatial profiles against test functions $\psi_i(z)$ to form weak matrix $\mathbf{G}$ and target $\mathbf{b}$.
2. **Evaluation & Constraint Grid:** Evaluates feature library terms across observational nodes and builds linear inequality system $A_{\text{ineq}} \boldsymbol{\xi} \ge \mathbf{b}_{\text{ineq}}$.
3. **Sparse Regression:** Solves for active coefficients using `STRidge` or `ConstrainedQP`.
4. **Intermediate Representation:** Wraps non-zero coefficients into `OperatorTerm{Float64}` objects and packages them inside a `DiscoveredModel` struct.

#### 2. Symbolic Closure Discovery (`discover_closure`)

Constructs a physically consistent `WSINDyClosure` directly from candidate `Vector{Num}` expressions:

```julia
result = discover_closure(
    data,
    candidates;
    target = :u,
    basis = GegenbauerBasis(n_spatial = 8, n_temporal = 1, lambda = 0.75),
    lambda = 1e-3,
    positivity_constraints = true
)

```

---

### Constrained Optimization & Closure Extraction

#### 1. Positivity-Constrained Sparse QP (`fit_wsindy_jump`)

When fitting closures, non-physical negative diffusivities ($K_m < 0$) are precluded by appending positivity bounds $\xi_j \ge 0$ for specified candidate indices $j \in \mathcal{I}_{\text{diff}}$:

$$\min_{\boldsymbol{\xi}} \, \frac{1}{2} \Vert{}\mathbf{G}\boldsymbol{\xi} - \mathbf{b}\Vert{}_2^2 + \lambda \Vert{}\boldsymbol{\xi}\Vert{}_1 \quad \text{subject to} \quad \mathbf{A}_{\text{ineq}} \boldsymbol{\xi} \ge \mathbf{b}_{\text{ineq}}$$

#### 2. Physical Closure Mapping (`extract_closure`)

Translates active candidate coefficients into usable closures for boundary layer integration:

$$K_m(z) = \sum_{k} \xi_k f_k(z)$$

$$K_h(z) = \frac{K_m(z)}{\operatorname{Pr}_t}, \qquad \operatorname{Pr}_t = 0.74$$

$$\tau_{\text{flux}} = u_*^2$$

---

### Component Summary

| Component | Function / Struct | Description |
| --- | --- | --- |
| **`GegenbauerBasis`** | `struct` | Spatial mode count $N_s$, temporal order $N_t$, and Gegenbauer weight parameter $\lambda$. |
| **`discover`** | `Function` | Primary API entry point mapping `ObservationTable` and `FeatureLibrary` to `DiscoveredModel`. |
| **`discover_closure`** | `Function` | Specialized pipeline returning weak matrices, coefficients, `DiscoveredModel`, and `WSINDyClosure`. |
| **`build_weak_library`** | `Function` | Projects candidate expressions against Gegenbauer test bases to generate $(\mathbf{G}, \mathbf{b})$. |
| **`fit_wsindy_jump`** | `Function` | Assembles inequality constraints ($A \boldsymbol{\xi} \ge d$) and calls `ConstrainedQP` solver via JuMP/HiGHS. |
| **`extract_closure`** | `Function` | Maps candidate expressions and sparse coefficient vectors into a functional `WSINDyClosure`. |