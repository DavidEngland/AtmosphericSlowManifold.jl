`src/Discretization/Backends/SpectralBLGalerkin.jl` provides a high-order Galerkin spectral discretization backend using Gegenbauer orthogonal polynomial expansions $C_n^{(\lambda)}(x)$. It transforms continuous 1D atmospheric boundary layer PDEs into low-dimensional mass-matrix ODE systems suitable for fast-slow manifold reduction, modal budget analysis, and stiff time integration.

---

### Mathematical Foundation & Spectral Projection

#### 1. Spatial Domain Mapping & Basis Polynomials

The physical height $z \in [0, H]$ is mapped onto the reference interval $x \in [-1, 1]$ via:


$$x(z) = \frac{2z}{H} - 1, \qquad z(x) = \frac{H}{2}(x + 1)$$

Spectral state profiles $u(t, z)$ are expanded in Gegenbauer modes $C_n^{(\lambda)}(x)$ under orthogonal weight $w_\lambda(x) = (1 - x^2)^{\lambda - 1/2}$:


$$u(t, z) \approx \sum_{n=1}^{N_{\text{modes}}} a_n(t) C_{n-1}^{(\lambda)}(x(z))$$

#### 2. Weak-Form Operators

* **Mass Matrix ($M$):**

$$M_{i,j} = \int_{-1}^1 C_{i-1}^{(\lambda)}(x) \, C_{j-1}^{(\lambda)}(x) \, (1 - x^2)^{\lambda - 0.5} \, dx$$


* **Stiffness Matrix ($K$):**

$$K_{i,j} = \int_0^H \frac{d C_{i-1}}{dz} \, \frac{d C_{j-1}}{dz} \, w_\lambda(x(z)) \, dz$$



#### 3. Mass-Matrix ODE System

The weak-form projection yields a stiff ODE system with mass matrix $M$:


$$M \frac{d\mathbf{a}}{dt} = -K \mathbf{a} + M \, \mathbf{N}(\mathbf{a})$$


where $\mathbf{N}(\mathbf{a})$ captures modal advection and diffusion interaction fluxes.

---

### Key Data Structures & System Components

| Component | Type / Signature | Description |
| --- | --- | --- |
| **`SpectralBLGalerkin`** | `AbstractDiscretization` | Primary discretization configuration ($N_{\text{modes}}$, weight parameter $\lambda$, domain height $H$, nonlinear scale switches). |
| **`SpectralNonlinearTensors`** | `struct` | Precomputed 3D Gegenbauer triple-product arrays (`triple`, `advection`, `diffusion_flux`) for zero-allocation modal nonlinear interactions. |
| **`BoundaryLayerWorkspace`** | `struct` | Real-space grid buffers (`z_grid`, `K_m_buffer`, `K_h_buffer`, `dK_dz_buffer`) coupling turbulence closures to spectral space. |
| **`ModalBudgetDiagnostic`** | `struct` | Decomposed modal tendencies ($\mathbf{f}_{\text{lin}}, \mathbf{f}_{\text{adv}}, \mathbf{f}_{\text{diff}}, \mathbf{f}_{\text{tot}}$) for transport budget tracking. |

---

### Non-Linear Interactions & Precomputed Tensors

To eliminate real-space quadrature inside inner ODE evaluation loops, Galerkin projections of non-linear products are precomputed into 3D tensors over $N_{\text{quad}}$ quadrature points:

$$\text{triple}_{k,i,j} = \int_{-1}^1 C_k C_i C_j \, w_\lambda(x) \, dz$$

$$\text{advection}_{k,i,j} = \int_{-1}^1 C_k C_i \frac{d C_j}{dz} \, w_\lambda(x) \, dz$$

$$\text{diffusion\_flux}_{k,i,j} = \int_{-1}^1 -\frac{d C_k}{dz} C_i \frac{d C_j}{dz} \, w_\lambda(x) \, dz$$

---

### Numerical Time Integration & Analytical Jacobians

When `dispatch_solve` is invoked for `SpectralBLGalerkin`:

1. **Closure Workspace Initialization:** Pre-allocates grid workspace and evaluates diffusivity scale parameters ($\text{adv\_scale}, \text{diff\_scale}$).
2. **Mass-Matrix Assembly:** Builds symmetric matrices $M$ and $K$ using numerical quadrature.
3. **Exact Analytical Jacobian (`jac_mass!`):** Derives exact Jacobian $\mathbf{J} = \frac{\partial \mathbf{f}}{\partial \mathbf{a}} = -K + M \frac{\partial \mathbf{N}}{\partial \mathbf{a}}$ for stiff implicit solvers.
4. **Stiff Solver Default:** Solves the resulting differential system using `Rodas5P()`, an $L$-stable 5th-order stiff Rosenbrock-Wanner scheme.

---

`src/Discretization/Backends/SpectralBLGalerkin.jl` provides a high-order Galerkin spectral discretization backend using Gegenbauer orthogonal polynomial expansions $C_n^{(\lambda)}(x)$. It transforms continuous 1D atmospheric boundary layer PDEs into low-dimensional mass-matrix ODE systems suitable for fast-slow manifold reduction, modal budget analysis, and stiff time integration.

---

### Mathematical Foundation & Spectral Projection

#### 1. Spatial Domain Mapping & Basis Polynomials

The physical height $z \in [0, H]$ is mapped onto the reference interval $x \in [-1, 1]$ via:

$$x(z) = \frac{2z}{H} - 1, \qquad z(x) = \frac{H}{2}(x + 1)$$

Spectral state profiles $u(t, z)$ are expanded in Gegenbauer modes $C_n^{(\lambda)}(x)$ under orthogonal weight $w_\lambda(x) = (1 - x^2)^{\lambda - 1/2}$:

$$u(t, z) \approx \sum_{n=1}^{N_{\text{modes}}} a_n(t) C_{n-1}^{(\lambda)}(x(z))$$

#### 2. Weak-Form Operators

* **Mass Matrix ($M$):**

$$M_{i,j} = \int_{-1}^1 C_{i-1}^{(\lambda)}(x) \, C_{j-1}^{(\lambda)}(x) \, (1 - x^2)^{\lambda - 0.5} \, dx$$

* **Closure-Weighted Stiffness Matrix ($K$):**

$$K_{i,j} = \int_0^H K_m(z) \, \frac{d C_{i-1}}{dz} \, \frac{d C_{j-1}}{dz} \, w_\lambda(x(z)) \, dz$$

Where the vertical eddy diffusivity profile $K_m(z)$ is interpolated via 1D linear interpolation (`_interp1d`) from closure buffer evaluation nodes directly onto the quadrature points $z(x)$. If no profile is provided, $K_m(z) \equiv 1.0$.

#### 3. Mass-Matrix ODE System

The weak-form projection yields a stiff ODE system with mass matrix $M$:

$$M \frac{d\mathbf{a}}{dt} = -K \mathbf{a} + M \, \mathbf{N}(\mathbf{a})$$

where $\mathbf{N}(\mathbf{a})$ captures modal advection and diffusion interaction fluxes.

---

### Key Data Structures & System Components

| Component | Type / Signature | Description |
| --- | --- | --- |
| **`SpectralBLGalerkin`** | `AbstractDiscretization` | Primary discretization configuration ($N_{\text{modes}}$, weight parameter $\lambda$, domain height $H$, nonlinear scale switches). |
| **`SpectralNonlinearTensors`** | `struct` | Precomputed 3D Gegenbauer triple-product arrays (`triple`, `advection`, `diffusion_flux`) for zero-allocation modal nonlinear interactions. |
| **`BoundaryLayerWorkspace`** | `struct` | Real-space grid buffers (`z_grid`, `K_m_buffer`, `K_h_buffer`, `dK_dz_buffer`) coupling physical $K_m(z)$ profiles directly into the spectral operator. |
| **`ModalBudgetDiagnostic`** | `struct` | Decomposed modal tendencies ($\mathbf{f}_{\text{lin}}, \mathbf{f}_{\text{adv}}, \mathbf{f}_{\text{diff}}, \mathbf{f}_{\text{tot}}$) for transport budget tracking. |

---

### Non-Linear Interactions & Precomputed Tensors

To eliminate real-space quadrature inside inner ODE evaluation loops, Galerkin projections of non-linear products are precomputed into 3D tensors over $N_{\text{quad}}$ quadrature points:

$$\text{triple}_{k,i,j} = \int_{-1}^1 C_k C_i C_j \, w_\lambda(x) \, dz$$

$$\text{advection}_{k,i,j} = \int_{-1}^1 C_k C_i \frac{d C_j}{dz} \, w_\lambda(x) \, dz$$

$$\text{diffusion\_flux}_{k,i,j} = \int_{-1}^1 -\frac{d C_k}{dz} C_i \frac{d C_j}{dz} \, w_\lambda(x) \, dz$$

---

### Numerical Time Integration & Analytical Jacobians

When `dispatch_solve` is invoked for `SpectralBLGalerkin`:

1. **Closure Workspace Initialization:** Pre-allocates grid workspace, updates $K_m(z)$ buffers via `update_diffusivity_buffers!`, and evaluates mean scalar response scales ($\text{adv\_scale}, \text{diff\_scale}$).
2. **Spatially Varying Matrix Assembly:** Computes mass matrix $M$ and evaluates stiffness matrix $K(K_m(z))$ by mapping real-space closure diffusivity values $K_m(z)$ onto quadrature nodes.
3. **Exact Analytical Jacobian (`jac_mass!`):** Derives exact mass-matrix Jacobian $\mathbf{J} = \frac{\partial \mathbf{f}}{\partial \mathbf{a}} = -K + M \frac{\partial \mathbf{N}}{\partial \mathbf{a}}$ for stiff implicit solvers.
4. **Stiff Solver Default:** Solves the resulting differential system using `Rodas5P()`, an $L$-stable 5th-order stiff Rosenbrock-Wanner scheme with autodiff Jacobians.