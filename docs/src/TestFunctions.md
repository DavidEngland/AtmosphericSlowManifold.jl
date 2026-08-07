`src/Discovery/TestFunctions.jl` implements smooth, compact-support test function families $\phi(z, t)$ for the Weak-form Sparse Identification of Nonlinear Dynamics (WSINDy) pipeline. Integrating PDE state variables against these test functions transfers derivative operators onto $\phi(z, t)$ via integration by parts, bypassing direct numerical differentiation of noisy field data.

---

### Test Function Families & Dispatches

| Struct / Function | Mathematical / Structural Purpose | Key Features & Boundary Conditions |
| --- | --- | --- |
| **`GegenbauerFamily`** | Global orthogonal polynomial test function family parameterized by $\lambda$ and `max_mode`. | Multiplied by mollifier $(1 - \xi^2)^2$ to enforce vanishing boundary conditions $\phi(z_0) = \phi(H) = 0$ at domain edges. |
| **`BSplineFamily`** | Localized compactly supported kernel family defined by polynomial `order` and grid `num_knots`. | Evaluates local radial bump functions $(1 - r)^{\text{order}}$ centered on knot positions, vanishing outside local kernel width. |
| **`evaluate_test_function`** | Computes spatial value $\phi(z)$ on domain $z \in [z_0, H]$. | Maps physical coordinate $z$ to normalized reference space $\xi \in [-1, 1]$. |
| **`evaluate_dt_test_function`** | Computes temporal derivative $\partial_t \phi(t)$ over integration window $t \in [t_0, T]$. | Polynomial time-scaling derivative using non-dimensional time $\tau = \frac{t - t_0}{T - t_0}$. |
| **`evaluate_dz2_test_function`** | Computes second spatial derivative $\partial_{zz} \phi(z)$. | Evaluates $\partial_{zz} \phi(z)$ via central finite differences with coordinate clamping `clamp(z, z0, H)`. |

---

### Mathematical Formulations

#### 1. Reference Coordinate Mapping

For physical spatial domain $z \in [z_0, H]$, coordinates are mapped onto reference domain $\xi \in [-1, 1]$ via:

$$\xi(z) = 2 \left( \frac{z - z_0}{H - z_0} \right) - 1$$

#### 2. Gegenbauer Test Function with Boundary Mollification

Gegenbauer modes $C_n^{(\lambda)}(\xi)$ are multiplied by a quartic envelope $(1 - \xi^2)^2$ to ensure smooth vanishing at boundary contacts $\xi = \pm 1$:

$$\phi_n(z) = C_n^{(\lambda)}(\xi(z)) \cdot \left( 1 - \xi(z)^2 \right)^2$$

#### 3. B-Spline Localized Support

For a given knot index $k$ with spacing $w$, radial coordinate $r$ determines local kernel strength:

$$r(z) = \frac{\vert \xi(z) - \xi_k \vert}{w}, \qquad \phi_k(z) = \begin{cases} (1 - r(z))^{\text{order}} & \text{if } r(z) < 1 \\ 0 & \text{otherwise} \end{cases}$$

#### 4. Second Spatial Derivative Numerical Approximation

Spatial second derivatives $\partial_{zz} \phi(z)$ are calculated using centered finite differences with step $h = 10^{-6}(H - z_0)$:

$$\frac{\partial^2 \phi}{\partial z^2}(z) \approx \frac{\phi(z + h) - 2\phi(z) + \phi(z - h)}{h^2}$$