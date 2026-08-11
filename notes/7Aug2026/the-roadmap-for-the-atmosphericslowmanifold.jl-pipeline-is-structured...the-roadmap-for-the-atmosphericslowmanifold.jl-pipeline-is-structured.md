The roadmap for the **AtmosphericSlowManifold.jl** pipeline is structured as a tiered progression that transitions from foundational data ingestion and symbolic discovery to advanced geometric analysis and probabilistic calibration [1-3]. The pipeline itself is built on a principle of **orthogonality**, decoupling observation space, manifold coordinates, symbolic closures, prognostic PDEs, and numerical discretization [4, 5].  
  
### 1. The Core Scientific Pipeline (Data Flow)  
The functional pipeline follows a logical hierarchy to ensure every discovered model remains physically interpretable [6].  
*   **Observation Layer:** Ingests heterogeneous field data (CASES-99, FLOSS, BLLAST, SHEBA) and automatically resolves aliases for friction velocity (\\(u_*\\)) and sensible heat flux (\\(H\\)) to derive the **Obukhov length scale (\\(L_{obukhov}\\))** [7-9].  
*   **Manifold & Geometry Engine:** Transforms raw physical variables into **intrinsic manifold coordinates** (modal amplitude \\(R\\), phase \\(\Omega\\), curvature \\(\chi\\), etc.) and analyzes invariant sets like fold curves and normal hyperbolicity [4, 10, 11].  
*   **Discovery Engine (WSINDy):** Projects governing PDEs onto smooth test function spaces to bypass noise, using **Sequential Thresholded Ridge Regression (STRidge)** to identify sparse symbolic operators [12-14].  
*   **System & Discretization:** Assembles the prognostic equations in **ModelingToolkit.jl** and solves them using either a **Stretched-Grid Finite Difference** backend or a **Modal Spectral Gegenbauer** backend featuring 3-tensor projections [15-17].  
  
### 2. Development Roadmap Tiers  
The development of the framework is categorized into four distinct tiers of implementation [1-3].  
  
| Tier | Status | Focus Area | Key Deliverables |  
| :--- | :--- | :--- | :--- |  
| **Tier 1** | **Completed** | Core Symbolic Engine | WSINDy discovery, weak-form assembly, and basic FD backends [1, 18]. |  
| **Tier 2** | **Completed** | Nonlinear Spectral Logic | 3-tensor Gegenbauer projections (\\(C_{ijk}^{(\lambda)}\\)), modal RHS operator decomposition, and solver-level smoke tests [19-21]. |  
| **Tier 3** | **Active** | GSPT & Bifurcation | Integration with **BifurcationKit.jl** to trace fold curves and Hopf points directly on discovered operators [2, 11, 22]. |  
| **Tier 4** | **Planned** | Advanced Calibration | **Hierarchical Bayesian Fitting** via Turing.jl to isolate site-specific parameters from global closures, and interactive human-in-the-loop dashboards [3, 23, 24]. |  
  
### 3. Immediate Functional Milestones  
To finalize the current active phase, the following actions are prioritized:  
*   **SHEBA Grid Reconstruction:** Resolving the vertical height interpolation for the SHEBA campaign to eliminate the current "N/A" transversality (\\(\mathcal{T}\\)) metric [3, 25].  
*   **Stiff Solver Hardening:** Optimizing the 12-hour prognostic benchmark using adaptive stiff solvers (e.g., **Rodas5P** or **RadauIIA5**) and applying \\(\tanh\\)-based stability saturation to prevent gradient blow-ups [16, 26, 27].  
*   **Automated Manuscript Generation:** Finalizing the **`make manuscript`** target, which automates the pipeline from raw data ingestion to the production of publication-ready LaTeX tables and figures [28-30].  
  
### 4. Long-Term Vision  
The final stage of the roadmap envisions the package as a **general-purpose discovery engine** for fast–slow geophysical systems [31]. This includes expanding into **Discrete Exterior Calculus (DEC)** or p-FEM backends for complex-terrain domains and developing a dedicated geometry subsystem for **canard detection** and desingularized flow computation [3, 32, 33].  
