`src/Geometry/CriticalManifold.jl` constructs the slow invariant set (critical manifold) $\mathcal{S}_0$ where fast turbulent transients have decayed and fast subsystem dynamics vanish ($\mathbf{g}(\mathbf{x}, \mathbf{y}) = 0$).

$$\mathcal{S}_0 = \left\{ (\mathbf{x}, \mathbf{y}) \in \mathbb{R}^m \times \mathbb{R}^n \;\middle\vert{}\; \mathbf{g}(\mathbf{x}, \mathbf{y}, \varepsilon = 0) = 0 \right\}$$

---

### Key Data Structures & Methods

#### Data Structures

* **`ManifoldPoint`**: Represents an individual point on $\mathcal{S}_0$, encapsulating slow variables ($\mathbf{x}$), fast variables ($\mathbf{y}$), and residual norm $\Vert{}\mathbf{g}(\mathbf{x}, \mathbf{y})\Vert{}$.
* **`CriticalManifoldSurface`**: Subtype of `AbstractInvariantSet` aggregating points across a slow-variable slice. Stores a coordinate array ($\mathbb{R}^{(m+n) \times N}$), fast nullclines, dimension metadata, and a `stability_mask` vector.

#### Algorithms

* **`find_manifold_point`**: Solves $\mathbf{g}(\mathbf{x}, \mathbf{y}) = 0$ for $\mathbf{y}$ given a fixed slow state $\mathbf{x}$ using a multi-variable Newton-Raphson iteration:

$$J = D_{\mathbf{y}}\mathbf{g}(\mathbf{x}, \mathbf{y}^k)$$


$$\mathbf{y}^{k+1} = \mathbf{y}^k - J^{-1}\mathbf{g}(\mathbf{x}, \mathbf{y}^k)$$


* **`solve_critical_surface`**: Traces the critical surface across a trajectory path $\mathbf{x}_{\text{path}}$ using warm-start seeding ($\mathbf{y}^0_{k+1} = \mathbf{y}_{k}$) to ensure fast convergence along continuous branches.

---

### Stability & Hyperbolicity Masking

During surface tracking in `solve_critical_surface`, the fast Jacobian spectrum is evaluated at each point:

```julia
J = finite_difference_jacobian_y(f_fast, p.x, p.y; params = params)
eigs = eigvals(Matrix{Float64}(J))
push!(stability_mask, all(real.(eigs) .< 0.0))

```

* **Attracting Branch (`stability_mask[i] == true`):** All eigenvalues satisfy $\operatorname{Re}(\lambda_j(D_{\mathbf{y}}\mathbf{g})) < 0$. The critical manifold is normally hyperbolic and attracting under Fenichel's Theorem.
* **Loss of Hyperbolicity / Fold Point (`stability_mask[i] == false`):** At least one eigenvalue satisfies $\operatorname{Re}(\lambda_j) \ge 0$. Transitions in the mask flag fold singularities ($\det D_{\mathbf{y}}\mathbf{g} = 0$), marking the physical boundary of turbulence collapse or regime change.

---

### Pipeline Integration

```
  Slow Path Iteration (x_path)
              │
              ▼
    find_manifold_point() ────> Newton-Raphson convergence to g(x, y) = 0
              │
              ▼
      Jacobian Spectrum ──────> Spectrum evaluation: Re(λ_j) < 0 ?
              │
              ▼
   CriticalManifoldSurface ───> Assembled invariant manifold + stability mask

```

This module provides the geometric base layer for `DesingularizedFlow.jl` (which rescales time along fold lines) and `FoldTracking.jl` (which traces the singular curves separating stable and unstable manifold sheets).