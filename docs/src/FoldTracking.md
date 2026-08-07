`FoldTracking.jl` provides the algorithms to track 1D and 2D fold lines ($\mathcal{F}$) where normal hyperbolicity breaks down ($\det D_{\mathbf{u}}\mathbf{f} = 0$). By formulating fold tracking as an augmented root-finding problem, it extracts both the singular locus, its tangent vectors, and local transversality metrics across varying background atmospheric parameters.

---

## Mathematical Formulation

A fold point $(\mathbf{u}^*, p^*)$ on a fast-slow manifold satisfies the augmented system condition $\mathbf{H}(\mathbf{u}, p) = \mathbf{0}$:

$$\mathbf{H}(\mathbf{u}, p) = \begin{bmatrix} \mathbf{f}(\mathbf{u}, p) \\ \det\left(D_{\mathbf{u}}\mathbf{f}(\mathbf{u}, p)\right) \end{bmatrix} = \mathbf{0}$$

Where:

* $\mathbf{f}(\mathbf{u}, p) = \mathbf{0}$ enforces that state $\mathbf{u}$ lies on the critical manifold $\mathcal{S}_0$.
* $\det\left(D_{\mathbf{u}}\mathbf{f}(\mathbf{u}, p)\right) = 0$ isolates the boundary where normal attraction vanishes.

### Tangent Vector & Transversality Condition

1. **Tangent Vector Calculation:** At the fold singularity, $D_{\mathbf{u}}\mathbf{f}$ has a zero eigenvalue. The nullspace direction $\mathbf{v} \in \ker(D_{\mathbf{u}}\mathbf{f})$ is isolated via Singular Value Decomposition (SVD):

$$D_{\mathbf{u}}\mathbf{f} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T \implies \mathbf{t} = \mathbf{V}[:, n]$$


2. **Transversality Metric ($\mathcal{T}$):** Measures the rate at which the determinant crosses zero with respect to parameter variation $p$:

$$\mathcal{T} = \frac{\partial}{\partial p} \det\left(D_{\mathbf{u}}\mathbf{f}(\mathbf{u}, p)\right)$$



Non-zero transversality ($\mathcal{T} \neq 0$) confirms a generic, structurally stable fold bifurcation rather than a higher-order degenerate singularity.

---

## Data Structure & Core API

### `FoldCurve` Struct

An `AbstractInvariantSet` subtype storing the continuation output along a parameter path:

| Field | Type | Description |
| --- | --- | --- |
| `points` | `Matrix{Float64}` | Array of size $(n+1) \times N$ containing concatenated state-parameter vectors $[\mathbf{u}; p]$. |
| `tangents` | `Matrix{Float64}` | Matrix of unit singular vectors spaning the zero-eigenvalue kernel at each fold point. |
| `transversality` | `Vector{Float64}` | Numeric derivative $\frac{d}{dp} \det(D_{\mathbf{u}}\mathbf{f})$ tracking crossing transversality. |

---

## Function Implementation Breakdown

### 1. `fold_indicator`

```julia
fold_indicator(f_fast, x, y; params)

```

Evaluates the local scalar determinant $\det(D_{\mathbf{y}}\mathbf{f})$ using central finite differences on $D_{\mathbf{y}}\mathbf{f}$. A value near zero flags proximity to a manifold fold boundary.

### 2. `_newton_scalar_correct`

```julia
_newton_scalar_correct(residual_fn, u0, p; tol, maxiter)

```

Provides a multidimensional Newton-Raphson solver tailored for augmented systems. It computes local Jacobians $J_{\mathbf{H}} = \frac{\partial \mathbf{H}}{\partial \mathbf{u}}$ via finite differences to correct candidate points back onto the constrained manifold:

$$\mathbf{u}^{k+1} = \mathbf{u}^k - J_{\mathbf{H}}^{-1} \mathbf{H}(\mathbf{u}^k, p)$$

### 3. `track_fold_line`

```julia
track_fold_line(residual_fn, u_init, p_range; n_steps, tol)

```

Executes parameterized fold continuation over a range $p \in [p_{\min}, p_{\max}]$:

1. Formulates the augmented residual vector $\mathbf{H}(\mathbf{u}, p) = [\mathbf{r}(\mathbf{u}, p); \det(J)]$.
2. Applies Newton correction to enforce $\mathbf{H}(\mathbf{u}, p) = \mathbf{0}$.
3. Computes the kernel vector via SVD of $J$ (`S.V[:, end]`).
4. Approximates transversality $\tau = \frac{d}{dp}\det(J)$ using a central difference derivative:

$$\tau \approx \frac{\det\left(J(p + h)\right) - \det\left(J(p - h)\right)}{2h}$$



### 4. `track_fold_curve`

```julia
track_fold_curve(f_fast, x_path, y_seed; params, tol)

```

Generates a `CriticalManifoldSurface` along a spatial path $\mathbf{x}_{\text{path}}$ using `solve_critical_surface`, filtering and returning all `ManifoldPoint` instances where $\vert{}\text{fold\_indicator}\vert{} \le \text{tol}$.

---

## Pipeline Role in Boundary Layer Dynamics

```
     Parameterized System r(u, p)
                 │
                 ▼
         track_fold_line()
                 │
  ┌──────────────┼──────────────┐
  ▼              ▼              ▼
Augmented     SVD Nullspace  Central Diff
Newton        Kernel         Transversality
  │              │              │
  └──────────────┼──────────────┘
                 │
                 ▼
         FoldCurve Output
 (Points, Tangents, Transversality)

```

In boundary layer transitions, `FoldTracking.jl` identifies parameter regimes (e.g., critical geostrophic wind speeds or surface cooling rates) where static stability collapses turbulent mixing modes, establishing the boundary between fully turbulent and strongly stable nocturnal layers.