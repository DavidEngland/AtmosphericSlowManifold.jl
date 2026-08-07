`src/Discovery/SymbolicExtraction.jl` defines the symbolic representation layer for model discovery in `AtmosphericSlowManifold.jl`. It converts numerical sparse regression outputs into `Symbolics.jl` and `ModelingToolkit.jl` (`Num`) expressions for symbolic verification, export, and simulation.

---

### Data Types & Structural Taxonomy

| Type | Supertype / Details | Description |
| --- | --- | --- |
| **`AbstractBasisFeature`** | Abstract Base | Root type for all candidate basis features in the library. |
| **`StateVariable`** | `<: AbstractBasisFeature` | Primary state variable (e.g., `:u`, `:v`, `:theta`). |
| **`SpatialDerivative`** | `<: AbstractBasisFeature` | Spatial derivative operator of order $k$ (e.g., variable `:u`, order `1` $\to \frac{\partial u}{\partial z}$). |
| **`DiagnosticVariable`** | `<: AbstractBasisFeature` | Auxiliary diagnostic field (e.g., diffusivity `:Km`, Richardson number `:Ri`). |
| **`BasisOperator`** | `struct` | Pairs a feature with a power exponent $p$ (`power::Float64`), representing $f^p$. |
| **`OperatorTerm{T}`** | `struct` | Represents a product term $c \prod_{j} f_j^{p_j}$ with coefficient $c::T$ and basis list `Vector{BasisOperator}`. |
| **`DiscoveredModel{T}`** | `struct` | Complete discovered PDE right-hand side equation with `target_variable`, term list, `residual_norm`, and `sparsity_level`. |

---

### Symbolic Mapping & MTK Expression Assembly

1. **Derivative Naming Convention:**
`get_feature_symbolic` maps spatial derivatives into symbolic lookup keys formatted as `d<order>_<variable>` (e.g., `SpatialDerivative(:u, 2)` maps to `:d2_u` in `var_map`).
2. **Monomial Reconstruction:**
`to_mtk_expression(term, var_map)` builds exact `Symbolics.Num` terms:

$$\text{Term} = c \times \prod_{b \in \text{basis}} b.\text{feature}^{b.\text{power}}$$


3. **Full PDE Model Assembly:**
`to_mtk_expression(model, var_map)` sums all active `OperatorTerm` instances into a single `Symbolics.Num` expression representing the discovered right-hand-side differential operator $\mathcal{N}(u, v, \theta, \dots)$.