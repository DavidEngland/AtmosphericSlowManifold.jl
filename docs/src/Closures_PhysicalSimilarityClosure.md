<!-- Auto-generated from package source -->
> **Source:** `src/Closures/PhysicalSimilarityClosure.md`

`src/Closures/PhysicalSimilarityClosure.jl` implements a data-driven similarity closure parameterized directly by observational campaign data or polynomial fits discovered via similarity diagnostics. It bridges empirical data collection (JSON campaign exports) with prognostic PDE backends while maintaining $C^\infty$ regularity.

---

### Core Mathematical & Algorithmic Design

1. **Horner's Method Polynomial Evaluation:**

$$\phi_m(\zeta) = \sum_{i=0}^N a_i \zeta^i$$



`_poly_eval` uses FMA (`muladd`) instructions in descending order for numerical stability and zero temporal overhead:
```julia
acc = muladd(acc, x, coeffs[i])

```


2. **Asymptotically Bounded Stability Metric ($\zeta$):**

$$\zeta = 5 \tanh\left( \frac{P_\zeta(z / L)}{5} \right)$$


* Prevents polynomial divergence ($\zeta \to \pm\infty$) at high vertical grid nodes ($z \gg L$) or under extreme stability regimes.
* `tanh` saturation bounds $\zeta \in [-5, 5]$ while remaining strictly $C^\infty$ smooth.


3. **$C^\infty$ Regularity Clamping:**

$$\phi_m(\zeta) = \operatorname{smooth\_floor}\left(P_\phi(\zeta), \; 0.1; \; \epsilon = 10^{-3}\right)$$


* Floor clamping at $0.1$ prevents singular diffusivities ($K_m \to \infty$) without introducing derivative step discontinuities.
* Height clamping at $z = 0$ uses `smooth_floor(z, 0)` for smooth surface layer bounds.


4. **Campaign JSON Serialization Ingestion:**
The `PhysicalSimilarityClosure(json_path::String)` constructor reads campaign diagnostic artifacts (`*_model_and_diagnostics.json`), parsing mean statistics for friction velocity ($u_*$), Obukhov length ($L$), and empirical polynomial terms with fallback resolution strategies.

---

### API Method Summary

| Function | Output / Expression | Mathematical Context |
| --- | --- | --- |
| `_km(c, z)` | $\frac{\kappa \cdot u_* \cdot z}{\phi_m(\zeta)}$ | Local momentum diffusivity ($m^2/s$). |
| `(closure)(z)` | $K_m(z)$ | Direct functor call interface for height $z$. |
| `eddy_momentum(c, m)` | $K_m(z_{\text{eval}})$ | Evaluates momentum diffusivity for a given `ManifoldState`. |
| `eddy_heat(c, m)` | $K_m(z_{\text{eval}}) / 1.1$ | Evaluates scalar heat diffusivity using a default Prandtl surrogate ($\mathrm{Pr}_t \approx 1.1$). |
| `surface_flux(c, m)` | $u_*^2$ | Kinematic surface momentum stress ($\tau_s$). |
| `evaluate_diffusivity_profile!` | $\mathbf{K}_m \in \mathbb{R}^{N_z}$ | In-place, zero-allocation profile evaluator for vertical momentum diffusion profiles. |
| `evaluate_heat_diffusivity_profile!` | $\mathbf{K}_h \in \mathbb{R}^{N_z}$ | In-place, zero-allocation profile evaluator for vertical heat diffusion profiles. |