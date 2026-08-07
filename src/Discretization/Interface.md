`src/Discretization/Interface.jl` defines the abstract dispatch interface for solving Single Column Model (SCM) boundary layer PDEs, establishing a strict decoupling between physical PDE definitions, turbulence closures, and spatial/temporal discretization strategies.

---

### Interface Architecture & Dispatch Pattern

1. **Facaded API (`solve_scm`):**
`solve_scm` serves as the primary top-level entry point for executing numerical integrations. It accepts symbolic `ModelingToolkit.PDESystem` definitions, an `AbstractClosure` representation, a discretization strategy subtype `disc::AbstractDiscretization`, and integration time span $t \in [t_0, t_f]$.
2. **Multiple Dispatch Hook (`dispatch_solve`):**
`dispatch_solve` provides the extensible backend protocol. Concrete subtypes of `AbstractDiscretization` must specialize `dispatch_solve` to convert symbolic MTK systems into numerical ODE/DAE problems and execute time integration via `DifferentialEquations.jl`.
3. **Fallback Error Handling:**
The default `dispatch_solve` method throws a clear `MethodError` if an unsupported or non-implemented discretization strategy is passed into `solve_scm`.

---

### Concrete Implementation Pattern

To implement a concrete discretization driver (such as a `MethodOfLines.jl`-backed spatial discretizer), specialize `AbstractDiscretization` and extend `dispatch_solve`:

```julia
# Example concrete implementation using MethodOfLines.jl
struct MOLDiscretization <: AbstractDiscretization
    z_nodes::Int
    order::Int
end

function MOLDiscretization(; z_nodes::Int = 100, order::Int = 2)
    return MOLDiscretization(z_nodes, order)
end

function dispatch_solve(
    disc::MOLDiscretization,
    pde_sys::PDESystem,
    closure::AbstractClosure,
    tspan::Tuple{Float64, Float64};
    alg = KenCarp4(),
    kwargs...
)
    # 1. Construct spatial grid specification
    grid_spec = [pde_sys.ivs[2] => disc.z_nodes] # ivs[2] corresponding to spatial coordinate 'z'
    discretization = MethodOfLines.MOLify(grid_spec; order = disc.order)

    # 2. Convert PDESystem to ODEProblem
    prob = MethodOfLines.discretize(pde_sys, discretization)

    # 3. Solve ODE using implicit stiff solver (e.g., KenCarp4, TRBDF2)
    sol = DifferentialEquations.solve(prob, alg; kwargs...)
    return sol
end

```

---

### Integration with Atmospheric Boundary Layer Pipeline

```
[ PrognosticPDE.jl ] ──> PDESystem ──┐
[ Closures.jl ]      ──> Closure    ──┼──> solve_scm(...) ──> dispatch_solve(...) ──> ODESolution
[ Interface.jl ]     ──> Strategy   ──┘

```

* **Stiffness Handling:** Implicit solvers (like `KenCarp4` or `TRBDF2`) leverage exact analytical Jacobians $\mathbf{J} = \frac{\partial \mathbf{f}}{\partial \mathbf{y}}$ derived symbolically from `ModelingToolkit.jl` through $C^\infty$-smooth closures.
* **Return Type:** Integration returns a standard `SciMLBase.ODESolution` containing spatial-temporal state profiles $u(z_i, t_j)$, $v(z_i, t_j)$, and $\theta(z_i, t_j)$.

Would you like to build out concrete discretization strategies (such as `MOLDiscretization` or a custom finite difference grid driver) in `src/Discretization/` next?

---

`src/Discretization/Interface.jl` is verified and fully aligned with the dispatch architecture.

### Interface Architecture

* **Abstract Base Type:** `AbstractDiscretization` provides the parent type for spatial discretization backends.
* **Dispatch Extension Point:** `dispatch_solve` enforces signature implementation across concrete discretization strategies, throwing a descriptive `MethodError` by default.
* **Public API Entry Point:** `solve_scm` acts as the standard pipeline call, forwarding all integration keyword arguments (`alg`, `reltol`, `abstol`, `saveat`) directly to backend solvers.

Would you like to implement a concrete spatial discretization driver (such as a `MethodOfLines.jl` or finite difference grid driver) in `src/Discretization/` next?