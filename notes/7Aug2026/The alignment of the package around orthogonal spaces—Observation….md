The alignment of the package around orthogonal spaces—**Observation**, **Manifold**, **Symbolic Closure**, **PDE System**, and **Discretization Backend**—provides a durable foundation. Moving closure discovery out of the numerical grid layer converts AtmosphericSlowManifold.jl from a standard Single-Column Model (SCM) into a general-purpose scientific discovery engine for fast–slow geophysical continuum systems.  
  
The 4-tier reprioritization cleanly separates the core identification machinery from downstream post-processing, and introducing a dedicated src/Geometry/ module formalizes GSPT as a first-class subsystem.  
  
## src/Geometry/ Subsystem Architecture  
Moving GSPT algorithms out of src/Manifold/ (which remains dedicated to state representations) and into src/Geometry/isolates geometric root-finding, manifold continuation, and hyperbolicity bounds into modular components.  
  
src/Geometry/  
├── Geometry.jl             # Submodule root & re-exports  
├── CriticalManifold.jl     # Root solver for 0 = f(x, y, 0) fast nullclines  
├── FoldTracking.jl         # Bordered system solver for det(D_y f) = 0 along parameter paths  
├── DesingularizedFlow.jl   # Rescaled slow flow & desingularized vector field on S_0  
├── CanardDetection.jl      # Folded singularity classification (nodes, saddles, foci)  
└── Fenichel.jl             # Normal hyperbolicity persistence threshold metrics  
**Primary Component Contracts**  
1. **CriticalManifold.jl** Computes the 3D critical manifold surface $S_0$ where fast boundary-layer dynamics reach quasi-equilibrium:  $$0 = f(\mathbf{x}, \mathbf{y}, 0)$$  Exposes non-linear root-finding interfaces (ManifoldPoint, CriticalManifoldSurface) parameterized by manifold coordinates $(R, \chi, \Pi_G)$.   
2. **FoldTracking.jl** Tracks fold curves $L \subset S_0$ satisfying $\det(D_{\mathbf{y}} f(\mathbf{x}, \mathbf{y}, 0)) = 0$. Uses ModelingToolkit.jl analytical Jacobians to construct bordered systems:  $$\begin{bmatrix} D_{\mathbf{y}} f & \mathbf{v} \\ \mathbf{w}^T & 0 \end{bmatrix} \begin{bmatrix} \mathbf{z} \\ \sigma \end{bmatrix} = \begin{bmatrix} \mathbf{0} \\ 1 \end{bmatrix}$$  yielding exact fold locations where $\sigma = 0$.   
3. **DesingularizedFlow.jl** Computes the slow flow on $S_0$:  $$D_{\mathbf{y}} f(\mathbf{x}, \mathbf{y}, 0) \cdot \dot{\mathbf{y}} = g(\mathbf{x}, \mathbf{y}, 0)$$  Applies the time-rescaling $d\tau = \det(D_{\mathbf{y}} f)^{-1} dt$ to generate the smooth desingularized vector field:  $$\dot{\mathbf{y}} = \operatorname{adj}(D_{\mathbf{y}} f) \cdot g(\mathbf{x}, \mathbf{y}, 0)$$   
4. **CanardDetection.jl** Identifies and classifies folded singularities where the desingularized vector field vanishes on the fold line $L$, separating trajectory trajectories that pass smoothly through fold points (canards).   
5. **Fenichel.jl** Evaluates normal hyperbolicity metrics by computing real parts of the fast Jacobian spectrum $\operatorname{Re}(\lambda_i(D_{\mathbf{y}} f))$ to establish spectral gaps and persistence bounds $\epsilon_0$.   
## Tier 1 Execution Plan: src/Discovery/WSINDyEngine.jl  
WSINDyEngine.jl provides sparse operator discovery over weak space-time formulations.  
  
[ Raw / Filtered Tower & LES Data ]  
                │  
                ▼  
  [ Space-Time Test Function Engine ]  <-- Evaluates test function integrals (G, b)  
                │  
                ▼  
    [ JuMP STRidge & Convex QP ]       <-- Min ||G Ξ - b||₂ + λ||Ξ||₁  s.t. K_m, K_h ≥ 0  
                │  
                ▼  
     [ Discovered WSINDyClosure ]      <-- Re-injected into build_pde_system()  
**Module Mechanics**  
1. **Space-Time Integration Engine** Avoids direct numerical differentiation of noisy boundary-layer profile data by integrating against smooth test functions $\phi_k(z, t) = \psi_i(z) \omega_j(t)$:  $$\int_{t_1}^{t_2} \int_{z_0}^{H} \phi_k \frac{\partial u}{\partial t} \, dz \, dt = - \int_{t_1}^{t_2} \int_{z_0}^{H} \frac{\partial \phi_k}{\partial t} u \, dz \, dt$$  $$\int_{t_1}^{t_2} \int_{z_0}^{H} \phi_k \frac{\partial}{\partial z} \left( K_m \frac{\partial u}{\partial z} \right) dz \, dt = \int_{t_1}^{t_2} \int_{z_0}^{H} \frac{\partial^2 \phi_k}{\partial z^2} K_m u \, dz \, dt$$   
2. **JuMP Constrained Optimization** Formulates Sequential Thresholded Ridge Regression (STRidge) and quadratic programs in JuMP.jl enforcing physical positivity constraints ($K_m(\mathbf{z}) \ge 0, K_h(\mathbf{z}) \ge 0$):  $$\min_{\mathbf{\Xi}} \frac{1}{2} \Vert{}\mathbf{G}\mathbf{\Xi} - \mathbf{b}\Vert{}_2^2 + \lambda \Vert{}\mathbf{\Xi}\Vert{}_1 \quad \text{s.t.} \quad \mathbf{A} \mathbf{\Xi} \ge \mathbf{0}$$   
3. **Symbolic AST Extraction** Translates non-zero sparse coefficients $\mathbf{\Xi}$ back into ModelingToolkit expressions, returning a fully runnable WSINDyClosure.   
## Tier 2 & 3 Roadmap Integration  
* **Tier 2 (Nonlinear Spectral Projections):** Upgrade src/Discretization/Backends/SpectralBLGalerkin.jl to compute 3-tensor inner products $\langle G_k^{(\lambda)}, G_i^{(\lambda)} G_j^{(\lambda)} \rangle_{\lambda}$ for advective terms $\left( u \frac{\partial u}{\partial z} \right)$ and non-linear closures $K_m(u, \theta) \frac{\partial u}{\partial z}$.   
* **Tier 3 (Bifurcation & Continuation):** Connect ModelingToolkit.jl Jacobians directly to BifurcationKit.jl. Because build_pde_system outputs a standard MTK system, bifurcations in discovered closures (e.g., transitions from turbulent mixing to laminar decoupling under strong stability) can be traced automatically along parameter curves like geostrophic wind $U_g$ or surface cooling rate $\dot{\theta}_0$.   
Shall we begin implementing **Tier 1** (src/Discovery/WSINDyEngine.jl) alongside the dedicated **src/Geometry/** module structure?  
