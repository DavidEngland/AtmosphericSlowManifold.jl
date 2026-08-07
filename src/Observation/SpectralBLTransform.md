`project_to_gegenbauer` projects discrete vertical atmospheric observation profiles onto an orthogonal **Gegenbauer spectral polynomial basis** $C_n^{(\lambda)}(x)$, providing compact, noise-filtered modal representations of boundary layer variables ($u, v, \theta, q$) for slow manifold diagnostic pipeline analysis and feature library integration.

---

### Mathematical Foundations

1. **Domain Normalization ($z \to x \in [-1, 1]$):**

$$x(z) = 2 \cdot \frac{z - z_{\min}}{z_{\max} - z_{\min}} - 1$$



Maps physical observation heights $z \in [z_{\min}, z_{\max}]$ onto the standard orthogonal polynomial domain $x \in [-1, 1]$.
2. **Three-Term Recurrence Relation ($C_n^{(\lambda)}$):**
Gegenbauer (ultraspherical) polynomials are generated iteratively in $O(n)$ time via `_obs_gegenbauerC`:

$$C_0^{(\lambda)}(x) = 1, \qquad C_1^{(\lambda)}(x) = 2\lambda x$$


$$C_k^{(\lambda)}(x) = \frac{2(k + \lambda - 1)x C_{k-1}^{(\lambda)}(x) - (k + 2\lambda - 2)C_{k-2}^{(\lambda)}(x)}{k}$$


3. **Weighted Galerkin Inner Product:**
Modal expansion coefficients $a_n$ are obtained by solving the weighted orthogonal projection:

$$a_n = \frac{\langle y, C_{n-1}^{(\lambda)} \rangle_{\lambda}}{\langle C_{n-1}^{(\lambda)}, C_{n-1}^{(\lambda)} \rangle_{\lambda}}$$


$$\langle f, g \rangle_{\lambda} = \int_{-1}^1 f(x) g(x) w(x; \lambda) \, dx, \qquad w(x; \lambda) = (1 - x^2)^{\lambda - 1/2}$$



---

### Key Numerical Safeguards

* **Singularity-Free Weight Clamping:** At domain boundaries $x = \pm 1$, the weighting function $(1 - x^2)^{\lambda - 1/2}$ can approach $0$ or encounter floating-point underflow. `_trapz_inner` clamps the weight function at a lower bound:

$$w(x; \lambda) = \max\left(10^{-12}, \; (1 - x^2)^{\lambda - 1/2}\right)$$


* **Physical Parameter Parameterization ($\lambda = 0.75$):**
* $\lambda = 0.50$: Chebyshev polynomials of the first kind / Legendre polynomials.
* $\lambda = 1.00$: Chebyshev polynomials of the second kind ($U_n$).
* $\lambda = 0.75$: Balances surface layer gradient resolution near $x = -1$ ($z = z_{\min}$) with bulk boundary layer weighting in the interior.



---

### Output NamedTuple Structure

Calling `project_to_gegenbauer(obs; n_modes=12, lambda=0.75)` returns:

| Field | Type | Description |
| --- | --- | --- |
| `coefficients` | `Dict{Symbol, Vector{Float64}}` | Modal coefficient vectors $a \in \mathbb{R}^{N_{\text{modes}}}$ for present profile variables (`:u`, `:v`, `:theta`, `:q`). |
| `lambda` | `Float64` | Ultraspherical parameter $\lambda$ used for the projection. |
| `n_modes` | `Int` | Number of spectral expansion modes retained ($N_{\text{modes}}$). |
| `z_span` | `Tuple{Float64, Float64}` | Min and max physical height boundaries $(z_{\min}, z_{\max})$. |
| `variables` | `Vector{Symbol}` | List of variable keys successfully projected into the basis. |