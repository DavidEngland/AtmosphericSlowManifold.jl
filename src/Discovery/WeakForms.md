`src/Discovery/WeakForms.jl` constructs the linear system matrix $\mathbf{G}\boldsymbol{\xi} \approx \mathbf{b}$ for Weak-form SINDy (WSINDy). By projecting observed boundary layer trajectories and candidate basis features onto smooth test function spaces $\phi(z, t)$, it eliminates the need to differentiate noisy observational data directly.

---

### Mathematical Formulation & Weak Projection Protocol

#### 1. Test Function Space & Space-Time Tensor Basis

Space-time test functions $\phi_{i,j}(z, t)$ are constructed as tensor products of spatial test modes $\psi_i(z)$ (from `GegenbauerFamily` or `BSplineFamily`) and temporal basis functions $\omega_j(\tau)$:

$$\phi_{i,j}(z, t) = \psi_i(z) \cdot \tau(t)^{j-1}, \qquad \tau(t) = 2 \left( \frac{t - t_{\min}}{t_{\max} - t_{\min}} \right) - 1$$

where $i \in \{1, \dots, N_{\text{spatial}}\}$ and $j \in \{1, \dots, N_{\text{temporal}}\}$.

#### 2. Weak Projection System Assembly

For target observation vector $\mathbf{y}$ and candidate basis library $f_c(z, t)$, row $r = (i-1)N_{\text{temporal}} + j$ of target vector $\mathbf{b}$ and candidate matrix $\mathbf{G}$ are assembled via 1D spatial trapezoidal integration over domain $z \in [z_0, H]$:

$$b_r = \int_{z_0}^H \phi_{i,j}(z, t) \cdot y(z, t) \, dz$$

$$G_{r,c} = \int_{z_0}^H \phi_{i,j}(z, t) \cdot f_c(z, t) \, dz$$

$$\int_{z_0}^H f(z) \, dz \approx \sum_{k=1}^{N-1} \frac{f(z_k) + f(z_{k+1})}{2} (z_{k+1} - z_k)$$

---

### Data Structures & Key Functions

| Entity | Type / Signature | Functional Role |
| --- | --- | --- |
| **`WeakFormMatrix`** | `struct` | Output container storing assembled weak matrix `G`, target vector `b`, and feature map array `feature_map`. |
| **`_normalize_axis`** | `Vector{Float64} -> Vector{Float64}` | Maps non-degenerate time coordinates $t$ into canonical reference interval $\tau \in [-1, 1]$. |
| **`_trapz`** | `(x, y) -> Float64` | Computes 1D composite trapezoidal numerical quadrature over non-uniform spatial or temporal grids. |
| **`_finite_diff_col`** | `(col, z, order) -> Vector{Float64}` | Computes $k$-th order central finite differences along vertical coordinate grid $z$ when evaluating spatial derivative features. |
| **`_candidate_value`** | `(expr::Num, i, z, data) -> Float64` | Symbolically substitutes sample states and evaluates candidate `Symbolics.Num` expressions numerically. |

---

### Assembly Methods & Overload Dispatches

1. **`FeatureLibrary` Overload:**
Used for standard WSINDy pipeline runs. Evaluates candidate features directly from `ObservationTable` columns or finite-difference spatial derivative arrays prior to weak-form trapezoidal projection.
2. **`Vector{Num}` Overload:**
Used for symbolic validation and custom candidate expressions. Evaluates candidate `Symbolics.Num` expressions at each vertical sample $i$ before projecting into test function space.