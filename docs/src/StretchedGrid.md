`src/Discretization/StretchedGrid.jl` implements a hyperbolic tangent coordinate transformation to construct non-uniform 1D vertical grids $z \in [0, H]$ with fine spatial resolution clustered near the surface boundary ($z = 0$).

---

### Mathematical Transformation

Given a uniform computational coordinate $s \in [0, 1]$ discretized into $N$ nodes:

$$s_i = \frac{i - 1}{N - 1}, \qquad i \in \{1, 2, \dots, N\}$$

The physical vertical coordinates $z_i \in [0, H]$ are calculated via:

$$z(s) = H \cdot \frac{\tanh(\alpha s)}{\tanh(\alpha)}$$

#### Grid Spacing & Local Resolution ($\Delta z$)

The local grid stretching factor (cell size derivative) is:

$$\frac{dz}{ds} = H \cdot \frac{\alpha \operatorname{sech}^2(\alpha s)}{\tanh(\alpha)}$$

* **Near-Surface Layer ($s \to 0$):** $\operatorname{sech}^2(0) = 1$, yielding maximum node density and minimal spacing $\Delta z_{\text{min}} \approx \frac{H}{N-1} \frac{\alpha}{\tanh(\alpha)}$.
* **Free Atmosphere ($s \to 1$):** Spacing expands exponentially to $\Delta z_{\text{max}} \approx \frac{H}{N-1} \frac{\alpha \operatorname{sech}^2(\alpha)}{\tanh(\alpha)}$.

---

### Parameter Specifications

| Parameter | Type | Range / Condition | Physical & Numerical Role |
| --- | --- | --- | --- |
| **`N`** | `Int` | $N \ge 2$ | Total number of vertical grid nodes along the column. |
| **`H`** | `Float64` | $H > 0$ | Domain top height $z_{\text{top}}$ (meters). |
| **`alpha` ($\alpha$)** | `Float64` | $\alpha > 0$ | Non-dimensional stretching intensity parameter. |

---

### Key Properties

1. **Exact Boundary Alignment:**

$$z(0) = 0, \qquad z(1) = H$$


2. **Uniform Limit ($\alpha \to 0^+$):**
Using L'Hôpital's rule:

$$\lim_{\alpha \to 0^+} \frac{\tanh(\alpha s)}{\tanh(\alpha)} = s \implies z(s) = H s \quad \text{(Uniform Grid)}$$


3. **Boundary Layer Resolution:**
Concentrates grid points in the lowest $100\text{ m}$ to resolve sharp surface-layer gradients in horizontal wind velocity ($\partial_z u, \partial_z v$) and potential temperature ($\partial_z \theta$) without incurring the computational cost of a globally fine uniform grid.