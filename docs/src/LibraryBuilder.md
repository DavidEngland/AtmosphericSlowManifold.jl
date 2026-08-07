`src/Discovery/LibraryBuilder.jl` constructs candidate feature libraries for symbolic model discovery, mapping physical state variables, spatial derivatives, and turbulence diagnostics into functional evaluators for candidate matrix generation in WSINDy.

---

### Library Generation Architecture

| Function / Struct | Inputs / Fields | Operational Role |
| --- | --- | --- |
| **`FeatureLibrary`** | `features::Vector{AbstractBasisFeature}`<br>

<br>`evaluators::Vector{Function}` | Container storing basis feature definitions and corresponding sample evaluation closures. |
| **`_feature_key`** | `f::AbstractBasisFeature` | Dispatches lookup keys matching sample dictionary symbols (`:u`, `:theta`, `:d1_u`, `:d2_v`, `:Km`). |
| **`build_feature_library`** | `states::Vector{Symbol}`<br>

<br>`diagnostics::Vector{Symbol}`<br>

<br>`max_derivative_order::Int` | Generates state features, spatial derivatives up to order $k$, diagnostic features, and builds type-cast `Float64` evaluation closures. |

---

### Feature Enumeration & Symbol Key Mapping

For a given configuration with states $S = \{\text{u}, \text{theta}\}$, diagnostics $D = \{\text{Km}\}$, and derivative order $K = 2$, `build_feature_library` builds the following feature array:

1. **State & Derivatives ($s \in S$):**
* `StateVariable(:u)` $\to$ key `:u`
* `SpatialDerivative(:u, 1)` $\to$ key `:d1_u` ($\partial_z u$)
* `SpatialDerivative(:u, 2)` $\to$ key `:d2_u` ($\partial_{zz} u$)
* `StateVariable(:theta)` $\to$ key `:theta`
* `SpatialDerivative(:theta, 1)` $\to$ key `:d1_theta` ($\partial_z \theta$)
* `SpatialDerivative(:theta, 2)` $\to$ key `:d2_theta` ($\partial_{zz} \theta$)


2. **Diagnostics ($d \in D$):**
* `DiagnosticVariable(:Km)` $\to$ key `:Km`



Each feature evaluator closure accepts sample dictionaries (`Dict{Symbol, Any}`) and extracts `float(sample[key])`, throwing a `KeyError` if required spatial field data or derivatives are missing.