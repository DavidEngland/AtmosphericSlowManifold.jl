## AtmosphericSlowManifold.jl: Project Overview & Status  
AtmosphericSlowManifold.jl is a high-performance Julia framework for atmospheric boundary layer (ABL) slow-manifold identification, spectral PDE discretization (Gegenbauer–Galerkin), data-driven closure discovery (WSINDyClosure), and multi-campaign observational analytics across CASES-99, FLOSS, BLLAST, and SHEBA.  
  
## Key Subsystem Architecture  
AtmosphericSlowManifold.jl  
├── src/  
│   ├── Observation/        # Campaign profile ingestion, flux aliasing, derived L_obukhov  
│   ├── Closures/           # MOST, WSINDy, and physical similarity closure models  
│   ├── Discretization/     # Stretched grids & Gegenbauer-Galerkin spectral solvers  
│   ├── System/             # Prognostic PDE engines, boundary conditions, export utilities  
│   └── AtmosphericSlowManifold.jl  
├── scripts/  
│   ├── run_campaign_exports.jl        # Production multi-campaign export pipeline  
│   ├── validate_campaign_exports.jl   # Schema and artifact validation suite  
│   └── run_pde_closure_benchmark.jl   # [In Progress] PDE solver benchmark  
├── reports/generated/campaign_exports/ # Artifacts (CSV, JSON, NetCDF, PNG, TeX)  
└── test/                              # Unit, integration, and schema regression suite  
## Campaign Diagnostic & Model Discovery Matrix  
All four observation campaigns now route through surface-flux-aware ingestion with derived Obukhov scaling ($L_{\text{obukhov}}$), isolating physical similarity parameters ($\zeta = z/L$, $\phi_{\text{obs}}$) instead of abstract basis terms ($a_0 \dots a_3$).  
  

| Campaign | Profiles (N) | Identified Model Terms | Valid Obukhov Scaling (N) | Mean Stability (ζˉ ) | Mean Stability Function (ϕˉ obs ) | Status |
| -------- | ------------ | ----------------------------- | ------------------------- | -------------------- | --------------------------------- | --------- |
| CASES-99 | $6,538$ | $\\phi_{\\text{obs}}, \\zeta$ | $6,538$ | $0.2114$ | $1.8211$ | Validated |
| FLOSS | $70,796$ | $\\phi_{\\text{obs}}, \\zeta$ | $70,796$ | $0.3842$ | $2.5368$ | Validated |
| BLLAST | $5,600$ | $\\phi_{\\text{obs}}, \\zeta$ | $5,600$ | $0.1985$ | $1.7940$ | Validated |
| SHEBA | $2,273$ | $\\phi_{\\text{obs}}, \\zeta$ | $2,273$ | $0.3450$ | $1.5011$ | Validated |
  
Code snippet  
  
\begin{table}[htbp]  
  \centering  
  \caption{Multi-Campaign Physical Parameter Extraction and Diagnostics}  
  \label{tab:campaign-overview-status}  
  \begin{tabular}{lccccc}  
    \toprule  
    \textbf{Campaign} & \textbf{Observations ($N$)} & \textbf{Discovered Terms} & \textbf{Obukhov $N$} & \textbf{Mean $\bar{\zeta}$} & \textbf{Mean $\bar{\phi}_{\text{obs}}$} \\  
    \midrule  
    \textbf{CASES-99} & 6,538 & $\phi_{\text{obs}}, \zeta$ & 6,538 & 0.2114 & 1.8211 \\  
    \textbf{FLOSS} & 70,796 & $\phi_{\text{obs}}, \zeta$ & 70,796 & 0.3842 & 2.5368 \\  
    \textbf{BLLAST} & 5,600 & $\phi_{\text{obs}}, \zeta$ & 5,600 & 0.1985 & 1.7940 \\  
    \textbf{SHEBA} & 2,273 & $\phi_{\text{obs}}, \zeta$ & 2,273 & 0.3450 & 1.5011 \\  
    \bottomrule  
  \end{tabular}  
\end{table}  
## Pipeline & Tooling Integration Status  

| Command Target | Primary Operational Function | Output Deliverables | Build Status |
| ---------------------- | --------------------------------------------------------- | ------------------------- | ------------ |
| make campaign-export | Batch data extraction and diagnostic computation | CSV, NetCDF, JSON, PNGs | Pass |
| make campaign-summary | Renders report summary tables in terminal | campaign_summary.md, .tex | Pass |
| make campaign-validate | Enforces artifact structure, schema, and parameter counts | Validation report output | Pass |
| make test | Executes core unit and integration test suite | Julia test harness | Pass |
| make pde-benchmark | Evaluates discovered physical closures in PDE solver | pde_benchmark_summary.csv | Next Up |
  
****Active Development Roadmap****  
1. **PDE Closure Coupling (PhysicalSimilarityClosure):** Coupling discovered physical closures ($\phi_{\text{obs}}, \zeta$) into GegenbauerGalerkin spectral discretizations (SpectralBLGalerkin.jl) and PrognosticPDESystem workspace buffers for zero-allocation ODE time-stepping.   
2. **PDE Benchmark Executable (make pde-benchmark):** Creating scripts/run_pde_closure_benchmark.jl to simulate 12-hour boundary layer evolutions and quantify profile drift reduction compared to unclosed neutral baselines ($K_m = \kappa u_* z$).   
3. **Automated LaTeX Manuscript Assembly (make manuscript):** Wiring LaTeX table fragments and comparative figures into an automated manuscript compilation pipeline.  
