The distinction between formal mathematical proof obligations and applied computational verification provides a clear, actionable path toward publication in journals like *Geoscientific Model Development (GMD)*, *Journal of Atmospheric Sciences (JAS)*, or *Journal of Open Source Software (JOSS)*. Treating omitted formal derivations as scope limits rather than fundamental framework flaws correctly recalibrates the evaluation criteria toward numerical stability, empirical verification, and reproducible software design.  
Below is the integrated classification, technical expansion on missing mathematical issues, and a revised, prioritized development roadmap for **AtmosphericSlowManifold.jl (ASM.jl)**.  
## 1. Audit Finding Reclassification  

| Finding / Area | Audit Severity | Reclassified Scope | Publication Requirement |
| -------------------------------- | -------------- | ------------------------------ | ------------------------------------------------------------------------------------------------------------------------- |
| Synthetic WSINDy Recovery | High | Critical Methodological Gap | Mandatory for peer-reviewed release; must prove exact recovery of known PDEs under additive noise. |
| Cross-Campaign Transferability | High | Validation Gap | Mandatory; Leave-One-Site-Out Cross-Validation (LOSO-CV) prevents overfitting claims. |
| Uncertainty Quantification (UQ) | Medium | Reporting Defect | Report all scalar metrics as \\text{Median} \\pm \\text{IQR} or 95\\% credibility intervals. |
| Benchmark Hardware Context | Low | Documentation Defect | Specify CPU, Julia version, state dimension N, and allocation profiles in benchmark outputs. |
| Test Function Specifications | Medium | Documentation Defect | Detail support radius R, polynomial degree, overlap ratio, and quadrature rule. |
| Fast-Slow Separation (\\epsilon) | High | Empirical Justification Needed | Document physical timescale ratios (\\tau_{\\text{turb}} / \\tau_{\\text{synoptic}}) rather than proving singular limits. |
| Canard Existence Proofs | High | Out of Scope (Numerical Only) | State explicitly that candidate canards are detected via numerical desingularization. |
| Fenichel Theorem Proof | High | Out of Scope (Spectral Check) | Track real parts of fast Jacobian eigenvalues \\Re(\\lambda_i(\\mathbf{D}_{\\mathbf{y}}\\mathbf{g})) \\neq 0 empirically. |
| ForwardDiff Instability | High | Non-Issue / Docs | Verify C^\\infty algebraic regularization of denominators before AD evaluation. |
| Gegenbauer Basis Selection | Medium | Documentation Defect | Move existing mathematical rationale for boundary clustering and adjustable \\lambda into main text. |
  
**2. Unaddressed Mathematical & Physical Issues**  
Two critical mathematical areas missing from the initial audit must be addressed to ensure robust model discovery.  
**Structural Identifiability & Sparse Library Coherence**  
When candidate libraries \mathbf{\Theta}(\mathbf{x}, \mathbf{y}) contain collinear or highly correlated physics terms (e.g., functions of Ri, \zeta = z/L, and \phi_m), WSINDy optimization can yield non-unique coefficient vectors \boldsymbol{\xi}.  
* **Mutual Coherence Metric:** Compute the cross-correlation matrix of test-function-projected library features \mathbf{G} = \mathbf{\Phi\Theta}: \mu(\mathbf{G}) = \max_{i \neq j} \frac{\vert{}\mathbf{g}_i^T \mathbf{g}_j\vert{}}{\Vert{}\mathbf{g}_i\Vert{}_2 \Vert{}\mathbf{g}_j\Vert{}_2} High mutual coherence (\mu \to 1) indicates severe ill-conditioning, requiring library pruning prior to thresholded least-squares (STLSQ) or Lasso steps.  
* **SVD & Condition Number:** Report the condition number \kappa(\mathbf{G}) = \sigma_{\max}/\sigma_{\min} for candidate feature matrices.  
* **Parameter Covariance Matrix:** Compute posterior parameter covariance \mathbf{\Sigma}_{\boldsymbol{\xi}} = \sigma^2 (\mathbf{G}^T \mathbf{G})^{-1} to quantify structural trade-offs between discovered terms.  
**Physical Invariants & Realizability Constraints**  
Relying solely on positivity bounds (K_m \ge 0) is insufficient to prevent unphysical dynamics. Discovery routines should enforce:  
* **Galilean Invariance:** PDE terms must remain invariant under constant velocity frame shifts \mathbf{u} \to \mathbf{u} + \mathbf{U}_0.  
* **Dimensional Homogeneity:** Enforce strict unit constraints across candidate library columns before matrix projection.  
* **Reynolds Stress Tensor Realizability:** Enforce Schumann-Lumley realizability conditions on fast turbulent variables: \overline{u_i' u_j'} \succeq 0 \quad \text{(Positive Semi-Definite Matrix)}   
* **Energy Dissipation Monotonicity:** Enforce that viscous/turbulent dissipation operators strictly act as energy sinks: \int_{\mathcal{Z}} u \cdot \mathcal{D}(u) \, dz \le 0   
## 3. Prioritized Development Roadmap  
```
                    ┌─────────────────────────────────────────┐
                    │ Phase 1: Synthetic & Validation Baseline│
                    │  - Synthetic WSINDy recovery + noise    │
                    │  - Cross-campaign (LOSO-CV) validation  │
                    └────────────────────┬────────────────────┘
                                         │
                                         ▼
                    ┌─────────────────────────────────────────┐
                    │ Phase 2: Identifiability & Constraints  │
                    │  - Mutual coherence & condition checks  │
                    │  - Realizability & dimensional bounds   │
                    └────────────────────┬────────────────────┘
                                         │
                                         ▼
                    ┌─────────────────────────────────────────┐
                    │ Phase 3: Diagnostic & Reporting Rigor   │
                    │  - Empirical Jacobian spectral tracking │
                    │  - Hardware specs & UQ intervals in tables│
                    └─────────────────────────────────────────┘

```
**Phase 1: Controlled Synthetic Recovery & Transferability**  
1. **Synthetic Noise Stress Test:**  
    * Generate canonical boundary layer trajectories from a known analytical model (e.g., Ekman spiral with known K_m(z)).  
    * Corrupt signals with additive Gaussian noise \eta \in \{0\%, 1\%, 5\%, 10\%, 20\%\}.  
    * Measure True Positive Rate (TPR) of recovered terms and L_2 norm coefficient error \Vert{}\boldsymbol{\xi}_{\text{discovered}} - \boldsymbol{\xi}_{\text{true}}\Vert{}_2.  
2. **Leave-One-Site-Out Cross-Validation (LOSO-CV):**  
    * Train library weights \boldsymbol{\xi} on three campaigns (\{\text{CASES-99}, \text{FLOSS}, \text{SHEBA}\}) and evaluate predictive residual norms on the held-out fourth dataset (\text{BLLAST}).  
**Phase 2: Structural Identifiability & Physical Hard Constraints**  
1. Incorporate mutual coherence pruning (\mu(\mathbf{G}) < 0.85) into WSINDyEngine.jl.  
2. Incorporate mutual coherence pruning (\mu(\mathbf{G}) < 0.85) into WSINDyEngine.jl.  
3. Incorporate mutual coherence pruning (\mu(\mathbf{G}) < 0.85) into WSINDyEngine.jl.  
4. Formulate term selection via Constrained Quadratic Programming (CQP) to enforce dimensional homogeneity and energy dissipation inequalities strictly.  
**Phase 3: Empirical GSPT Checks & Diagnostic UQ**  
1. **Empirical Normal Hyperbolicity Tracking:**  
    * Compute and store minimum magnitude real parts of fast eigenvalues along integrated trajectories: \gamma(t) = \min_i \left\vert{} \Re\left(\lambda_i\left(\mathbf{D}_{\mathbf{y}}\mathbf{g}\right)\right) \right\vert{}   
    * Define empirical normal hyperbolicity loss as \gamma(t) < \delta_{\text{thresh}}.  
2. **Standardized UQ Reporting:**  
    * Replace single-point scalar metrics across all generated LaTeX tables with \text{Median} \,\, [Q_{0.025}, Q_{0.975}].  
3. **Benchmarking Environment Manifest:**  
    * Add automated hardware context logging (CPU microarchitecture, thread count, Julia/BLAS version, RAM allocations) directly to benchmark JSON/TeX exports.  
Which phase of this updated roadmap would you like to implement first into the repository suite?  
