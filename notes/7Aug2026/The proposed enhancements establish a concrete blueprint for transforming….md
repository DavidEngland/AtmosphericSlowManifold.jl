The proposed enhancements establish a concrete blueprint for transforming AtmosphericSlowManifold.jl into a production-grade scientific engine. Below is the technical specification and data architecture for implementing these five streams across the codebase.  
  
## 1. Advanced Geometric & Bifurcation Engine  
To bridge data-driven operator discovery with GSPT, the geometry subsystem (src/Geometry/) uses compiled Jacobians from ModelingToolkit.jl to drive pseudo-arclength continuation via BifurcationKit.jl.  
  
**Desingularized Vector Field & Canard Classification**  
Along the fold curve $\mathcal{F}(\mathbf{z}) = \{ \mathbf{z} \mid \det(D_{\mathbf{y}}\mathbf{g}(\mathbf{z})) = 0 \}$, the fast subsystem loses normal hyperbolicity. The desingularized vector field on the slow manifold is defined as:  
  
$$\frac{d\mathbf{x}}{d\tau} = \text{adj}\left(D_{\mathbf{y}}\mathbf{g}(\mathbf{y}, \mathbf{x})\right) \cdot \mathbf{f}(\mathbf{y}, \mathbf{x})$$  
A folded singularity occurs at points $\mathbf{z}^* \in \mathcal{F}(\mathbf{z})$ where the desingularized vector field vanishes:  
  
$$\text{adj}\left(D_{\mathbf{y}}\mathbf{g}(\mathbf{y}^*, \mathbf{x}^*)\right) \cdot \mathbf{f}(\mathbf{y}^*, \mathbf{x}^*) = \mathbf{0}$$  
Trajectories passing through $\mathbf{z}^*$ form canard solutions, crossing from stable to unstable manifold sheets without immediate turbulence collapse.  
  
**Output Data Formats**  
* **Bifurcation Branches:** reports/generated/bifurcations/branch_data.csv containing continuation parameter values, solution norms, leading eigenvalues $\lambda_{1,2}$, and bifurcation point flags (Fold, Hopf).   
* **Geometric Portraits:** Phase plots (.pdf/.png) displaying slow manifold geometry, fold curves, and canard passages.   
## 2. Hierarchical Bayesian Calibration & Uncertainty Quantification  
To distinguish universal similarity physics from site-specific boundary conditions, src/Closures/BayesianCalibration.jlformulates a multi-site hierarchical model in Turing.jl.  
  
**Probabilistic Hierarchy**  
$$\boldsymbol{\theta}_{\text{global}} \sim \text{MvNormal}(\boldsymbol{\mu}_0, \boldsymbol{\Sigma}_0) \quad \text{(Universal WSINDy coefficients)}$$  
$$z_{0, s} \sim \text{LogNormal}(\mu_{z0, s}, \sigma_{z0}^2) \quad \text{(Site-specific surface roughness for campaign } s \text{)}$$  
$$\phi_{\text{obs}, k, s} \sim \mathcal{N}\left( \phi_m\left(\zeta_{k, s}; \boldsymbol{\theta}_{\text{global}}\right), \sigma_{\text{model}}^2 + \sigma_{\text{site}, s}^2 \right)$$  
**Output Data Formats**  
* **MCMC Parameter Distributions:** reports/generated/calibration/posterior_samples.json holding posterior means, standard deviations, and credible intervals for cross-site transfer learning.   
* **Convergence Diagnostics:** ELBO trace logs, Effective Sample Size ($ESS$), and Gelman-Rubin $\hat{R}$statistics exported as LaTeX tables.   
## 3. Diagnostics & Computational Hygiene  
To ensure zero-allocation runtime performance and physical consistency, src/Diagnostics/ exposes dedicated tracking structures.  
  
**Modal Energy Budget Decomposition**  
The Gegenbauer modal kinetic energy equation for mode $k$ is tracked via:  
  
$$\frac{d E_k}{d t} = \mathcal{P}_k^{\text{shear}} + \mathcal{B}_k^{\text{buoyancy}} - \mathcal{D}_k^{\text{dissipation}}$$  
Julia  
  
struct ModalBudgetDiagnostic  
    timestamp::Float64  
    mode_index::Int  
    kinetic_energy::Float64  
    advective_transfer::Float64  
    diffusivity_dissipation::Float64  
    coriolis_work::Float64  
end  
**Output Data Formats**  
* **Budget Metrics:** reports/generated/diagnostics/modal_energy_budget.csv or .nc (NetCDF) containing modal power balances over time.   
* **Hygiene Logs:** Benchmark allocation logs confirming 0 Bytes allocated per RHS call in SpectralBLGalerkin.   
## 4. Interactive Interface & Automated Publication Pipeline  
**Interactive Discovery Dashboard**  
A web dashboard (built via Plotly/Dash in Julia) exposes:  
  
1. **Symbolic Term Pruning:** Real-time toggling of WSINDy candidate library terms with instant Pareto frontier update (Complexity vs. RMSE).   
2. **$Ri_{\text{cr}}$ Knee Visualizer:** Interactive phase diagrams of TKE versus flux Richardson number ($Ri$) with overlaid GSPT fold locations.   
**Reproducible Manuscript Pipeline (make manuscript)**  
The Makefile target orchestrates the execution chain from raw observational data to a compiled paper:  
  
Makefile  
  
.PHONY: manuscript  
  
manuscript:  
	@echo "==> Step 1/4: Ingesting Field Campaigns & Reconstructing Grids..."  
	julia --project=. -e 'using AtmosphericSlowManifold; DataIngestion.reprocess_all()'  
	@echo "==> Step 2/4: Running WSINDy Discovery & GSPT Geometry Engine..."  
	julia --project=. -e 'using AtmosphericSlowManifold; Discovery.run_pipeline()'  
	@echo "==> Step 3/4: Executing PDE Benchmark Suite..."  
	julia --project=. -e 'using AtmosphericSlowManifold; Benchmarks.run_pde_suite()'  
	@echo "==> Step 4/4: Generating LaTeX Artifacts and Compiling Paper..."  
	cd reports/paper && pdflatex main.tex  
## 5. SHEBA Vertical Grid Reconstruction  
To eliminate N/A metric entries and zero-height level artifacts, src/Observation/DataIngestion.jl interpolates SHEBA mast measurements onto standard vertical coordinates $z \in [2.0, 3.0, 5.0, 10.0, 18.0] \text{ m}$.  
  
Julia  
  
function standardize_sheba_profiles!(raw_df::DataFrame)  
    # Group observations by timestamp and interpolate height levels  
    standard_heights = [2.0, 3.0, 5.0, 10.0, 18.0]  
    grouped = groupby(raw_df, :timestamp)  
      
    processed_dfs = Vector{DataFrame}(undef, length(grouped))  
    for (i, group) in enumerate(grouped)  
        # Logarithmic interpolation for wind velocity and potential temperature  
        u_interp = LinearInterpolation(group.z, group.u, extrapolation_bc=Line())  
        theta_interp = LinearInterpolation(group.z, group.theta, extrapolation_bc=Line())  
          
        df_i = DataFrame(  
            timestamp = group.timestamp[1],  
            z = standard_heights,  
            u = u_interp.(standard_heights),  
            theta = theta_interp.(standard_heights)  
        )  
        processed_dfs[i] = df_i  
    end  
    return vcat(processed_dfs...)  
end  
This reconstruction allows the geometry engine to derive $L_{\text{obukhov}}$, compute $\zeta$, and calculate the transversality metric $\bar{\mathcal{T}}$ for SHEBA.  
  
Would you like to start by implementing the SHEBA grid reconstruction script to finalize the multi-campaign diagnostic summary table?  
