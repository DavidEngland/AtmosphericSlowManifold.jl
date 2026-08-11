## Discovery Subsystem Architecture (src/Discovery/)  
With the src/Geometry/ layer active, the next architectural milestone decomposes WSINDyEngine.jl into single-responsibility submodules centered around a structured **Intermediate Representation (IR)**.  
  
src/Discovery/  
├── Discovery.jl            # Submodule entrypoint & exports  
├── TestFunctions.jl        # Test function families (Gegenbauer, B-Spline, Fourier)  
├── WeakForms.jl            # Space-time quadrature & weak-form G, b matrix assembly  
├── LibraryBuilder.jl       # Candidate feature matrix construction (u, θ, Ri, Π_G, ∂u/∂z)  
├── ConstraintBuilder.jl    # Physical inequalities (K_m ≥ 0, monotonicity bounds)  
├── SparseRegression.jl     # Optimization backends (STRidge, ConstrainedQP, MIOSR)  
├── SymbolicExtraction.jl   # OperatorTerm IR & ModelingToolkit AST generation  
└── WSINDyEngine.jl         # Orchestrator & backward-compatible discover() entrypoint  
## Intermediate Representation: OperatorTerm & DiscoveredModel  
Before building ModelingToolkit ASTs, discovered coefficients are stored as structured metadata objects. This intermediate layer decouples sparse regression solvers from symbolic AST construction, enabling operator simplification, unit checking, and model export (JSON/YAML).  
  
Julia  
  
# Defined in src/Discovery/SymbolicExtraction.jl  
  
@enum OperatorKind DIFFERENTIAL ALGEBRAIC DIAGNOSTIC  
  
struct BasisOperator  
    symbol::Symbol             # e.g., :u, :theta, :Ri, :Pi_G, :dz_u  
    kind::OperatorKind  
    spatial_derivative_order::Int  
    power::Float64  
end  
  
struct OperatorTerm  
    coefficient::Float64  
    basis::Vector{BasisOperator}  
end  
  
struct DiscoveredModel  
    target_variable::Symbol  
    terms::Vector{OperatorTerm}  
    residual_norm::Float64  
    sparsity_level::Int  
end  
## Module Responsibilities & Execution Order  
1. **TestFunctions.jl**   
    * Implements AbstractTestFunctionFamily with GegenbauerFamily and BSplineFamily.   
    * Evaluates compact-support space-time test functions $\phi_k(z, t) = \psi_i(z) \omega_j(t)$ and their exact analytical derivatives $\frac{\partial \phi_k}{\partial t}$, $\frac{\partial^2 \phi_k}{\partial z^2}$.   
2. **WeakForms.jl**   
    * Performs 2D space-time integration-by-parts on profile data to populate weak linear system matrices $\mathbf{G}\mathbf{\Xi} \approx \mathbf{b}$:  $$\int_{t_1}^{t_2} \int_{z_0}^{H} \phi_k \frac{\partial u}{\partial t} \, dz \, dt = - \int_{t_1}^{t_2} \int_{z_0}^{H} \frac{\partial \phi_k}{\partial t} u \, dz \, dt$$   
3. **LibraryBuilder.jl**   
    * Evaluates candidate algebraic, differential, and diagnostic features across the spatio-temporal domain to populate candidate feature columns in $\mathbf{G}$.   
4. **ConstraintBuilder.jl**   
    * Evaluates candidate features over a dense evaluation grid to construct convex inequality matrices $\mathbf{A}_{\text{ineq}} \mathbf{\Xi} \ge \mathbf{b}_{\text{ineq}}$, enforcing physical bounds such as non-negative eddy diffusivities ($K_m \ge 0, K_h \ge 0$).   
5. **SparseRegression.jl**   
    * Provides optimization backends (STRidge and ConstrainedQP via JuMP.jl / HiGHS) to solve:  $$\min_{\mathbf{\Xi}} \frac{1}{2} \Vert{}\mathbf{G}\mathbf{\Xi} - \mathbf{b}\Vert{}_2^2 + \lambda \Vert{}\mathbf{\Xi}\Vert{}_1 \quad \text{s.t.} \quad \mathbf{A}_{\text{ineq}} \mathbf{\Xi} \ge \mathbf{b}_{\text{ineq}}$$   
6. **SymbolicExtraction.jl**   
    * Maps non-zero sparse coefficient vectors $\mathbf{\Xi}$ to OperatorTerm structures and translates them into executable ModelingToolkit.Num ASTs for WSINDyClosure.   
7. **WSINDyEngine.jl**   
    * Orchestrates the pipeline end-to-end while maintaining backward compatibility with existing discover_closurefunction signatures.   
## Planned Verification Tests  
* **test/test_wsindy_ir.jl**: Validates OperatorTerm construction, term simplification, and conversion into valid ModelingToolkit.Num expressions.   
* **test/test_wsindy_pipeline.jl**: Validates weak-matrix assembly, positivity constraint enforcement, and sparse recovery of synthetic turbulent viscosity profiles ($K_m(z) = k_0 + k_1 z$).   
Ready to implement the src/Discovery/ submodules and IR layer?  
