# AtmosphericSlowManifold.jl — Project Status Report  
AtmosphericSlowManifold.jl is scaffolded, precompiled, and passing all unit and integration tests. The package establishes a modular Julia framework that embeds Weak Sparse Identification of Non-Linear Dynamics (WSINDy) inside a Geometric Singular Perturbation Theory (GSPT) manifold engine, decoupling symbolic closure discovery from spatial discretization.  
  
## 1. Completed System Architecture  
                       [ Observation Space ]  
                Tower CSV & NetCDF Profile Readers  
                                 │  
                                 ▼  
                        [ Manifold Space ]  
           ManifoldState (η_i, R, Ω, χ, Π_G, λ_min)  
             & GSPT Fold Diagnostics (det(J_y) = 0)  
                                 │  
                                 ▼  
                     [ Symbolic Closure API ]  
           AbstractClosure: WSINDyClosure & MOSTClosure  
                                 │  
                                 ▼  
                    [ Prognostic PDE System ]  
             PrognosticPDE.jl (ModelingToolkit.jl)  
                                 │  
                 ┌───────────────┴───────────────┐  
                 ▼                               ▼  
       [ MethodOfLinesFD ]              [ SpectralBLGalerkin ]  
    Stretched Grid FD Backend        Modal ODE Gegenbauer Engine  
     (Operational / WRF SCM)              (ROM / GSPT Bifurcation)  
## 2. Core Capabilities Implemented  
* **Intrinsic Manifold State (src/Manifold/):**   
    * ManifoldState defines symbolic coordinates $(\eta_1, \eta_2, \eta_3, R, \Omega, \chi, \Pi_G, \lambda_{\min})$and physical state variables ($u, v, \theta, q, u_*, z_0$).   
    * GSPTDiagnostics.jl evaluates loss of normal hyperbolicity along slow manifold trajectories via fast Jacobian determinants $\det(J_{\mathbf{y}}) = 0$.   
* **Composable Symbolic Closures (src/Closures/):**   
    * AbstractClosure interface enforcing standard dispatch for eddy_momentum, eddy_heat, and surface_flux.   
    * WSINDyClosure holds discovered symbolic expressions; MOSTClosure provides Monin-Obukhov similarity baseline functions.   
* **Prognostic PDE Engine (src/System/):**   
    * PrognosticPDE.jl constructs ModelingToolkit.PDESystem conservation equations for momentum and potential temperature.   
    * Symbolic substitution maps manifold variables to state profile expressions before executing derivative expansions (expand_derivatives).   
* **Dual Discretization Backends (src/Discretization/):**   
    * MethodOfLinesFD: Stretched vertical grid generator using hyperbolic tangent clustering ($z_i = H \frac{\tanh(\alpha i / N)}{\tanh(\alpha)}$) wired to MethodOfLines.jl.   
    * SpectralBLGalerkin: Reduced-order Galerkin solver evaluating weighted Gegenbauer mass ($\mathbf{M}$) and stiffness ($\mathbf{K}$) matrices via recurrence relations.   
* **Strict Observation Ingestion (src/Observation/):**   
    * DataIngestion.jl enforces strict unit and column validation ($z\,[\text{m}], u, v, u_*\,[\text{m s}^{-1}], \theta\,[\text{K}], q\,[\text{kg kg}^{-1}]$) across CSV records and NetCDF profile files (radiosondes, Cabauw, NEON, LES outputs).   
## 3. Test & Verification Status  
Running julia --project=. -e 'using AtmosphericSlowManifold; include("test/runtests.jl")' yields green passes across all test modules:  
  
1. **test_wsindy_discovery.jl:** Confirms symbolic expression substitution, unit consistency, and operator dispatch for WSINDyClosure.   
2. **test_gspt_fold.jl:** Validates fold condition evaluation and hyperbolicity loss checks on manifold trajectories.   
3. **test_scm_backends.jl:** Proves **architecture invariance**—swapping closures (MOSTClosure $\leftrightarrow$WSINDyClosure) and backends (MethodOfLinesFD $\leftrightarrow$ SpectralBLGalerkin) executes without modifying the prognostic PDE signature or function calls.   
## 4. Immediate Roadmap & Pending Work  
1. **JuMP WSINDy Matrix Engine (src/Discovery/WSINDyEngine.jl):** Connect weak-form space-time integral matrices $\mathbf{G}\mathbf{\Xi} \approx \mathbf{b}$ to JuMP.jl for constrained STRidge sparse regression over Gegenbauer test function domains.   
2. **Full Nonlinear Spectral Projections (src/Discretization/Backends/SpectralBLGalerkin.jl):** Extend the Gegenbauer linear stiffness engine to project nonlinear advection and closure terms $\langle G_n^{(\lambda)}, \mathcal{N}(\mathbf{U}) \rangle_{\lambda}$.   
3. **Hierarchical Bayesian Site Calibration (src/Calibration/HierarchicalTuring.jl):** Complete the Turing.jl MCMC pipeline to separate global operator parameters $\mathbf{\Xi}_{\text{global}}$ from site-specific surface adapters ($z_0, z_{0h}, \text{LAI}, \Pi_G$).  
