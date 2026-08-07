`src/Discovery/Discovery.jl` acts as the primary module aggregator for the symbolic model discovery engine in `AtmosphericSlowManifold.jl`. It assembles the Weak-form Sparse Identification of Nonlinear Dynamics (WSINDy) framework, physical constraint enforcement, candidate library generation, and symbolic verification into a unified discovery pipeline.

---

### Discovery Submodule Architecture

| Inclusion File | Subsystem Functionality | Mathematical & Operational Role |
| --- | --- | --- |
| **`SymbolicExtraction.jl`** | Feature Extraction | Converts state profiles, spatial derivatives, and diagnostic closures into symbolic feature variables. |
| **`LibraryBuilder.jl`** | Candidate Library Construction | Builds candidate function matrices $\mathbf{\Theta}(\mathbf{U})$ containing polynomial combinations, spatial derivatives, and closure couplings. |
| **`ConstraintBuilder.jl`** | Physical Constraints | Enforces linear/equality boundary conditions, mass/momentum conservation rules, and coefficient bounds on candidate terms ($\mathbf{C} \boldsymbol{\xi} = \mathbf{d}$). |
| **`TestFunctions.jl`** | Compact Test Functions | Generates $C^\infty$ smooth test function families $\phi(z, t)$ (e.g., localized bump functions or Chebyshev mollifiers) with compact support. |
| **`WeakForms.jl`** | Weak Formulation Integration | Integrates space-time field data against test functions via integration-by-parts ($\int u \partial_t \phi \, dz dt$), eliminating noisy numerical differentiation. |
| **`SparseRegression.jl`** | Sparse Solvers | Implements Sequentially Thresholded Least Squares (STLS) and Sparse Relaxation (SR3) algorithms to select dominant terms. |
| **`WSINDyEngine.jl`** | Core Discovery Engine | Orchestrates the end-to-end weak SINDy pipeline, assembling linear system $\mathbf{G}\boldsymbol{\xi} \approx \mathbf{b}$ and solving for active terms. |
| **`SymbolicVerification.jl`** | Model Validation | Reconstructs symbolic MTK/Symbolics expressions from sparse coefficient vectors $\boldsymbol{\xi}$ and computes residual error norms. |

---

### WSINDy Model Discovery Pipeline

```
[ Raw Trajectories u(z,t) ] ──> [ TestFunctions.jl ] ──> [ WeakForms.jl ] ──┐
                                                                            ├──> [ WSINDyEngine.jl ] ──> [ SparseRegression.jl ] ──> DiscoveredModel
[ LibraryBuilder.jl ] ────────> [ Candidate Library Θ ] ─────────────────────┤                                    ▲
                                                                            │                                    │
[ ConstraintBuilder.jl ] ────> [ Equality Constraints Cξ = d ] ────────────┘                        [ ConstraintBuilder.jl ]

```

1. **Noise-Robust Weak Formulation:** Rather than differentiating noisy boundary layer observations directly, `WeakForms.jl` projects target field equations onto smooth test functions from `TestFunctions.jl`:

$$\int_{\Omega} \frac{\partial u}{\partial t} \phi_k \, dz dt = -\int_{\Omega} u \frac{\partial \phi_k}{\partial t} \, dz dt$$


2. **Constrained Sparse Regression:** `SparseRegression.jl` solves the penalized, constrained minimization problem over candidate library $\mathbf{G}$ and target vector $\mathbf{b}$:

$$\min_{\boldsymbol{\xi}} \Vert{}\mathbf{G}\boldsymbol{\xi} - \mathbf{b}\Vert{}_2^2 + \lambda \Vert{}\boldsymbol{\xi}\Vert{}_0 \quad \text{subject to} \quad \mathbf{C}\boldsymbol{\xi} = \mathbf{d}$$


3. **Symbolic Output Generation:** `SymbolicVerification.jl` maps non-zero coefficients back into readable algebraic expressions, returning a `DiscoveredModel` struct ready for export or simulation.