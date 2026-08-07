`src/Geometry/Jacobians.jl` provides symbolic differentiation, compiled kernel generation, and numeric fallback mechanisms for fast-subsystem Jacobians ($D_{\mathbf{y}}\mathbf{g}$), determinants, adjugate matrices, and singular tangent space projections.

---

### Core Architectural Mechanics

1. **Symbolic Compilation Pipeline**
* The constructor `JacobianModel(f_fast, fast_vars, slow_vars)` automatically computes symbolic Jacobians, determinants, and adjugate matrices via `Symbolics.jl`.
* `_compile_symbolic_kernel` uses `Symbolics.build_function(..., expression = Val(false))` to compile dynamic symbolic expressions directly into executable Julia functions (`J_function`, `det_function`, `adj_function`). This eliminates runtime symbolic evaluation overhead during high-frequency manifold tracking and numerical continuation.


2. **Analytical Adjugate Algebra for Desingularization**
* Explicit analytical cofactor matrices are implemented for $2 \times 2$ and $3 \times 3$ symbolic systems:

$$\text{adj}(J) = C^T, \quad C_{ij} = (-1)^{i+j} M_{ij}$$


* **GSPT Significance:** In Singular Perturbation Theory, inverse Jacobians $J^{-1} = \frac{\text{adj}(J)}{\det(J)}$ singularize at fold loci where $\det(D_{\mathbf{y}}\mathbf{g}) = 0$. Pre-computing $\text{adj}(D_{\mathbf{y}}\mathbf{g})$ enables the construction of desingularized slow vector fields:

$$\frac{d\mathbf{x}}{d\tau} = \text{adj}(D_{\mathbf{y}}\mathbf{g}) \cdot \mathbf{f}(\mathbf{x}, \mathbf{y})$$



allowing smooth numerical integration across singular manifold fold lines.


3. **Nullspace & Tangent Space Kernel Extraction**
* `evaluate_tangent_space(J)` uses Singular Value Decomposition (SVD) ($J = U \Sigma V^T$) to isolate right singular vectors corresponding to singular values below threshold ($s_i \le \text{thresh}$).
* The extracted columns of $V$ define the nullspace $\ker(D_{\mathbf{y}}\mathbf{g})$, supplying the exact kernel direction along which normal hyperbolicity breaks down.


4. **Dual Numeric-Symbolic Evaluation Strategy**
* **Compiled Model Path:** `evaluate_jacobian`, `evaluate_det`, and `evaluate_adjugate` execute in-memory compiled kernels on Float64 vectors (`u_fast`, `u_slow`, `p`).
* **Black-Box Finite Difference Path:** `compute_fast_jacobian(f_fast, x, y)` provides a 2nd-order central finite difference fallback for arbitrary Julia functions where symbolic derivatives are unavailable:

$$\frac{\partial f_i}{\partial y_j} \approx \frac{f_i(\mathbf{x}, \mathbf{y} + \varepsilon \mathbf{e}_j) - f_i(\mathbf{x}, \mathbf{y} - \varepsilon \mathbf{e}_j)}{2\varepsilon}$$





---

### API Reference & Method Summary

| Method | Inputs | Output | Mathematical / GSPT Function |
| --- | --- | --- | --- |
| `JacobianModel` | `f_fast`, `fast_vars`, `slow_vars`, `params` | `JacobianModel` struct | Builds and compiles symbolic $D_{\mathbf{y}}\mathbf{g}$, $\det(D_{\mathbf{y}}\mathbf{g})$, and $\text{adj}(D_{\mathbf{y}}\mathbf{g})$. |
| `compute_adjugate` | Matrix $J$ ($2 \times 2$, $3 \times 3$, or Numeric) | Transposed Cofactor Matrix | Analytical adjugate for symbolic systems; $\det(J) J^{-1}$ for numeric matrices. |
| `compute_fast_jacobian` | Model + `Dict` substitutions OR Function + numeric vectors | `Matrix{Float64}` | Evaluates Jacobian numerically via compiled kernels or central finite differences. |
| `evaluate_tangent_space` | Jacobian Matrix $J$, `atol` | `Matrix{Float64}` | Computes nullspace basis $\ker(J)$ via SVD for fold direction identification. |

---

Would you like to examine `CriticalManifold.jl` next to see how root-finding constructs the slow invariant manifold $\mathcal{S}_0$, or explore `DesingularizedFlow.jl` to see how the adjugate kernels are used for trajectory continuation through fold points?