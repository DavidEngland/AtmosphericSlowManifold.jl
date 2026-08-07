`ManifoldMetrics` provides geometric and dynamical systems metrics for characterizing slow-manifold stability, fold bifurcation proximity, fast-slow timescale separation, and off-manifold relaxation error in reduced atmospheric models.

---

### Geometric & Dynamical Metrics Summary

| Function | Mathematical Formulation | Return Range / Type | Geometric & Physical Interpretation |
| --- | --- | --- | --- |
| **`transversality`** | $$\tau = \frac{\vert \mathbf{v}_{\text{fast}} \cdot \mathbf{n}_{\text{slow}} \vert}{\Vert{}\mathbf{v}_{\text{fast}}\Vert{}_2 \Vert{}\mathbf{n}_{\text{slow}}\Vert{}_2}$$

 | $[0, 1]$ (`Float64`) | Measures normalized alignment between fast tendency vector $\mathbf{v}_{\text{fast}}$ and manifold normal vector $\mathbf{n}_{\text{slow}}$. Values near `1.0` indicate orthogonal, transverse fast approach to the slow manifold. |
| **`fold_distance`** | $$d(\mathbf{x}, \mathcal{F}) = \min_{j} \Vert{}\mathbf{x} - \mathbf{f}_j\Vert{}_2$$

 | $[0, \infty)$ (`Float64`) | Computes minimum Euclidean $L_2$ distance from state $\mathbf{x}$ to discrete points $\mathbf{f}_j$ along fold-locus matrix $\mathbf{F}$. Indicates proximity to fold bifurcations where manifold stability changes. |
| **`slow_manifold_error`** | $$\text{SME} = \Vert{}\mathbf{x} - \mathcal{M}(\mathbf{x})\Vert{}_2$$

 | $[0, \infty)$ (`Float64`) | Measures $L_2$ norm of spatial error between current state $\mathbf{x}$ and its projected state on slow-manifold mapping $\mathcal{M}(\mathbf{x})$. |
| **`normal_hyperbolicity`** | $$R = \frac{\vert\text{Re}(\lambda_{\text{fast}, n_{\text{fast}}})\vert}{\vert\text{Re}(\lambda_{\text{slow}, 1})\vert}$$

 | $[0, \infty)$ (`Float64`) | Computes spectral gap ratio using Jacobian eigenvalues sorted by descending magnitude of real components $\vert\text{Re}(\lambda)\vert$. Higher ratios ($R \gg 1$) signify strong timescale separation and manifold persistence. |

---

### Key Implementation Features

1. **Zero-Allocation Views in Distance Traversal:**
`fold_distance` uses `@views for i in eachindex(state)` to traverse matrix columns without creating heap allocations or array slices during high-frequency trajectory checks.
2. **Eigenvalue Spectral Ordering:**
`normal_hyperbolicity` extracts real components $\text{Re}(\lambda_i)$, sorts by descending magnitude (`sort(abs.(real.(vals)); rev = true)`), and computes spectral gap across boundary index $n_{\text{fast}}$.
3. **Division-by-Zero Safety:**
* `transversality` checks `iszero(den)` to handle degenerate zero-norm vectors, returning `0.0`.
* `normal_hyperbolicity` checks `iszero(slow_rate)` to prevent divide-by-zero errors when slow mode rates vanish, returning `typemax(Float64)`.
* `fold_distance` handles empty locus matrices ($m = 0$) by returning `typemax(Float64)`.