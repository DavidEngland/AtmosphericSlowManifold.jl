This code implements the classification and detection layer for **folded singularities** and **canard trajectories** in `AtmosphericSlowManifold.jl`.

In Geometric Singular Perturbation Theory (GSPT), folded singularities are points on the fold line ($\det D_{\mathbf{y}}\mathbf{g} = 0$) where the desingularized slow vector field vanishes. They dictate whether boundary layer state trajectories can cross smoothly from attracting (stable) to repelling (unstable) sheets of the slow invariant manifold without undergoing instantaneous turbulent breakdown jumps.

---

### Mathematical Classification Framework

The function `classify_folded_singularity` linearizes the desingularized vector field $\mathbf{w}(\mathbf{x}, \mathbf{y}) = \operatorname{adj}(D_{\mathbf{y}}\mathbf{g})\mathbf{f}(\mathbf{x},\mathbf{y})$ around a candidate fold point, evaluating the local Jacobian matrix $J = D_{\mathbf{y}}\mathbf{w}$. The local dynamics are categorized via the trace ($\operatorname{tr}(J)$), determinant ($\det(J)$), and discriminant ($\Delta = \operatorname{tr}(J)^2 - 4\det(J)$):

| Singular Type | Symbolic Enum | Algebraic Condition | Geometric & Canard Dynamics |
| --- | --- | --- | --- |
| **Folded Saddle** | `FOLDED_SADDLE` | $\det(J) < -tol$ | **True Canard:** Exactly two singular canard trajectories cross from the attracting sheet to the repelling sheet. |
| **Folded Node** | `FOLDED_NODE` | $\det(J) > tol, \; \Delta \ge 0$ | **Canard Funnel:** A continuous family (funnel) of trajectories passes smoothly across the fold line onto the repelling manifold. |
| **Folded Focus** | `FOLDED_FOCUS` | $\det(J) > tol, \; \Delta < 0$ | **No Canards:** Complex conjugate eigenvalues cause spiraling dynamics; trajectories cannot cross the fold line smoothly and undergo fast jumps. |
| **Degenerate** | `nothing` | $\vert{}\det(J)\vert{} \le tol$ | Higher-order or non-hyperbolic bifurcation requiring higher-order normal forms. |

---

### Detailed Component Analysis

#### 1. Data Structures (`FoldedSingularity` & `CanardSegment`)

* **`FoldedSingularity`**: Serves as a diagnostic container binding the base `ManifoldPoint`, string/symbol classification, enumerated `SingularType`, and exact complex eigenvalues of $J$.
* **`CanardSegment`**: An `AbstractInvariantSet` subtype that stores the concatenated state space coordinate $[x; y]$, the verified `SingularType`, the complex eigenvalue spectrum, and the numerically integrated trajectory matrix `canard_trajectory`.

#### 2. Vector Field Linearization (`detect_folded_singularity`)

Computes the local $n \times n$ Jacobian of the desingularized vector field (`vfield`) at a given `ManifoldPoint` using a second-order central finite difference scheme:

$$J_{ij} = \frac{\partial w_i}{\partial y_j} \approx \frac{w_i(\mathbf{x}, \mathbf{y} + \varepsilon \mathbf{e}_j) - w_i(\mathbf{x}, \mathbf{y} - \varepsilon \mathbf{e}_j)}{2\varepsilon}$$

#### 3. Segment Constructor (`build_canard_segment`)

Acts as a strict constructor guard. It asserts that `folded.singular_type !== nothing`, preventing degenerate or unclassified singularities from generating invalid invariant manifold segments.

---

### Atmospheric & Boundary Layer Significance

```
                   Slow Manifold S_0
      Attracting Sheet (Stable)    Repelling Sheet (Unstable)
     =========================\   /========================
                               \ /
                          Folded Saddle
                       (Canard Trajectory)
                                │
                                ▼
         Smooth, delayed transition across turbulence collapse
         without immediate fast-transient numerical instability.

```

In stable boundary layers (SBLs), the emergence of **folded saddles** and **folded nodes** explains persistent non-equilibrium states. Rather than immediately collapsing into laminar flow when passing a critical threshold (like $Ri_{\text{cr}}$), atmospheric trajectories following canard paths stay near the repelling manifold for $O(1)$ slow time scales. Capturing these trajectories allows `ASM.jl` to prevent spurious numerical oscillations during transient nocturnal transitions.