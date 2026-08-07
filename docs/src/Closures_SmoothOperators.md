<!-- Auto-generated from package source -->
> **Source:** `src/Closures/SmoothOperators.md`

`src/Closures/SmoothOperators.jl` provides generic, $C^\infty$-differentiable algebraic approximations for non-smooth threshold operators ($\max$, $\min$, lower-bound clamping). It eliminates non-differentiable $C^0$ kinks across closure models, allowing `ModelingToolkit.jl` and `Symbolics.jl` to construct continuous Jacobians and Hessians for continuation solvers and GSPT diagnostics.

---

### Mathematical Properties & Formulations

The operators use algebraic hyperbolic smoothing to approximate piecewise linear functions:

1. **Smooth Maximum ($\operatorname{smooth\_max}$):**

$$\operatorname{smooth\_max}(x; \epsilon) = \frac{x + \sqrt{x^2 + \epsilon^2}}{2} \approx \max(x, 0)$$


* Asymptotics: $\lim_{x \to \infty} \operatorname{smooth\_max}(x; \epsilon) = x$, $\lim_{x \to -\infty} \operatorname{smooth\_max}(x; \epsilon) = 0$, and $\operatorname{smooth\_max}(0; \epsilon) = \frac{\epsilon}{2}$.


2. **Smooth Minimum ($\operatorname{smooth\_min}$):**

$$\operatorname{smooth\_min}(x; \epsilon) = \frac{x - \sqrt{x^2 + \epsilon^2}}{2} \approx \min(x, 0)$$


3. **Smooth Floor Clamping ($\operatorname{smooth\_floor}$):**

$$\operatorname{smooth\_floor}(x, y; \epsilon) = y + \operatorname{smooth\_max}(x - y; \epsilon) \approx \max(x, y)$$



---

### Key Technical Characteristics

| Feature | Implementation Detail | Architectural Benefit |
| --- | --- | --- |
| **Type Promotion (`promote`)** | `promote(x, eps)` / `promote(x, floor_val, eps)` | Ensures type stability when mixing native floating-point types (`Float64`, `Float32`) with symbolic variables (`Symbolics.Num`). |
| **$C^\infty$ Regularity** | Smooth square-root algebraic form | Replaces conditional step transitions with infinitely differentiable curves, eliminating Dirac-delta derivative spikes. |
| **Zero-Allocation Inlining** | `@inline` compiler directives | Inlines calculations directly into spatial PDE discretization loops (`@simd`) with $0\text{ Bytes}$ memory allocation. |

---

### Integration Across the Solver Pipeline

```
          Raw Inputs (x, floor_val, eps)
                       │
                       ▼
          promote(x, floor_val, eps)
  (Handles Float64, Float32, Symbolics.Num)
                       │
                       ▼
            Hyperbolic Smoothing
        (f1 + (dx + sqrt(dx² + e²))/2)
                       │
         ┌─────────────┴─────────────┐
         ▼                           ▼
   PDE RHS Solvers            Symbolic ASTs
  (Zero allocation)      (ModelingToolkit.jl)
         │                           │
         ▼                           ▼
 Newton Continuations        Analytical Jacobians
 (Continuation.jl)          (GSPT / Fenichel Metrics)

```

1. **Newton Convergence:** Eliminates step discontinuities in the Jacobian matrix $D_{\mathbf{y}}\mathbf{f}$ that otherwise cause Newton-Raphson solvers in `Continuation.jl` (`continue_manifold`) to overshoot or fail near stability bounds.
2. **Fenichel Manifold Persistence:** Guarantees that fast-slow vector fields satisfy the $C^r$ smoothness ($r \ge 1$) required by Fenichel's Theorem, preventing spurious non-hyperbolic artifacts in spectral boundary calculations (`fenichel_metrics`).