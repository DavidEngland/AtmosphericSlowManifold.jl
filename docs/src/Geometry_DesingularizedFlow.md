<!-- Auto-generated from package source -->
> **Source:** `src/Geometry/DesingularizedFlow.md`

These functions implement the time-rescaled slow flow dynamics and folded singularity detection for `DesingularizedFlow.jl`. They resolve numerical blowups near singular fold manifolds ($\det D_{\mathbf{y}}\mathbf{f} = 0$) by applying matrix adjugate rescaling.

---

### Mathematical Formulations

#### 1. Standard Slow Flow (`slow_flow_vector`)

On the critical manifold $\mathcal{S}_0$, fast dynamics vanish ($\mathbf{f}(\mathbf{x}, \mathbf{y}) = 0$). Differentiating this constraint along slow trajectories yields the constrained slow flow:

$$D_{\mathbf{y}}\mathbf{f}(\mathbf{x}, \mathbf{y}) \frac{d\mathbf{y}}{dt} = -D_{\mathbf{x}}\mathbf{f}(\mathbf{x}, \mathbf{y}) \mathbf{g}(\mathbf{x}, \mathbf{y})$$

$$\frac{d\mathbf{y}}{dt} = -\left(D_{\mathbf{y}}\mathbf{f}\right)^{-1} D_{\mathbf{x}}\mathbf{f} \cdot \mathbf{g}(\mathbf{x}, \mathbf{y})$$

* **Singularity:** As trajectories approach the fold manifold where $\det(D_{\mathbf{y}}\mathbf{f}) = 0$, the standard slow flow vector diverges ($J \backslash g \to \infty$).

#### 2. Desingularized Flow (`desingularized_vector_field`)

To integrate through or analyze flow behavior at fold points, the independent variable is rescaled via $d\tau = \frac{dt}{\det(D_{\mathbf{y}}\mathbf{f})}$. Using Cramer’s rule ($J^{-1} = \frac{\text{adj}(J)}{\det(J)}$), multiplying by $\det(D_{\mathbf{y}}\mathbf{f})$ removes the singular denominator:

$$\frac{d\mathbf{y}}{d\tau} = \operatorname{adj}\left(D_{\mathbf{y}}\mathbf{f}\right) \mathbf{g}(\mathbf{x}, \mathbf{y})$$

This vector field remains smooth and bounded everywhere on $\mathcal{S}_0$, including on the fold line $\mathcal{F}$.

---

### Folded Singularity Identification

The function `find_desingularized_singular_points` scans candidate points on `CriticalManifoldSurface` to identify **folded singularities** (e.g., folded saddles, nodes, or foci):

$$\begin{cases} \det\left(D_{\mathbf{y}}\mathbf{f}(\mathbf{x}, \mathbf{y})\right) = 0 & \text{(Fold Locus Condition)} \\ \operatorname{adj}\left(D_{\mathbf{y}}\mathbf{f}(\mathbf{x}, \mathbf{y})\right) \mathbf{g}(\mathbf{x}, \mathbf{y}) = \mathbf{0} & \text{(Desingularized Equilibrium Condition)} \end{cases}$$

```
                Critical Manifold S_0
              ┌───────────────────────┐
              │  det(J) > 0 (Attracting)
              │                       │
 Fold Line ───┼───────────────────────┼──> det(J) = 0
              │                       │
              │  det(J) < 0 (Repelling)
              └───────────┬───────────┘
                          │
                          ▼
            [ Folded Singularity Point ]
        det(J) ≈ 0  AND  ||adj(J)·g|| ≈ 0
                          │
                          ▼
             Allows Canard Trajectories
            (Passage Stable -> Unstable)

```

### Key Numerical Properties

* **Threshold Checking:** A point is classified as a folded singularity if both the determinant indicator $d = \vert{}\det(J)\vert{} \le \text{det\_tol}$ and the rescaled velocity $\Vert{}\mathbf{v}\Vert{} = \Vert{}\operatorname{adj}(J)\mathbf{g}\Vert{} \le \text{flow\_tol}$ fall within user-defined tolerances.
* **Canard Precursor:** Isolating these points is a prerequisite for `CanardDetection.jl`, as true canard trajectories (which cross smoothly from attracting to repelling manifold sheets without undergoing fast jumps) can only pass through folded singularities.