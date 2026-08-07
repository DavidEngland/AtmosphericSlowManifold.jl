AtmosphericSlowManifold.jl: A Data-Driven Framework for Boundary Layer Dynamics

Executive Summary

AtmosphericSlowManifold.jl (ASM.jl) is a high-performance Julia framework designed to revolutionize the modeling and analysis of the atmospheric boundary layer (ABL). By unifying Geometric Singular Perturbation Theory (GSPT) with Weak-form Sparse Identification of Non-linear Dynamics (WSINDy), the framework replaces traditional, empirical stability functions—which often fail during rapid non-equilibrium transitions—with data-driven, physically constrained differential equations.

Critical Takeaways:

* Geometric Transition Criteria: Unlike classical models that rely on static thresholds like the critical Richardson number (Ri_{\text{cr}}), ASM.jl identifies turbulence collapse as a loss of normal hyperbolicity on a critical slow manifold.
* Noise-Robust Discovery: The WSINDy engine utilizes weak-form integration against C^\infty test functions, allowing for the symbolic discovery of governing equations without the noise amplification inherent in numerical differentiation.
* Computational Efficiency: The framework employs Spectral Galerkin discretization with precomputed 3-tensor contractions, enabling zero-allocation right-hand-side (RHS) evaluations essential for stiff ODE solvers.
* Hierarchical Bayesian Calibration: ASM.jl utilizes hierarchical pooling to estimate parameters across heterogeneous field campaigns, balancing universal physical constants with site-specific variations like surface roughness.
* Validated Performance: The system has been benchmarked against major field datasets (CASES-99, FLOSS, BLLAST, SHEBA), demonstrating the ability to capture complex dynamics like sunset boundary layer collapse and extreme polar stability.

1. Core Philosophy: From Empirical to Geometric

Classical turbulence closure models, specifically Monin–Obukhov Similarity Theory (MOST), rely on empirical stability functions \phi_m(\zeta) and \phi_h(\zeta). These assume local thermodynamic equilibrium—a condition that breaks down during:

* Sunset boundary layer collapse.
* Low-level jet formation.
* Extreme stability over polar ice sheets.

ASM.jl addresses these failures by modeling the ABL as a fast-slow dynamical system: \frac{d\mathbf{x}}{dt} = \mathbf{f}(\mathbf{x}, \mathbf{y}, \epsilon), \qquad \epsilon \frac{d\mathbf{y}}{dt} = \mathbf{g}(\mathbf{x}, \mathbf{y}, \epsilon)

* \mathbf{x} (Slow variables): Mean wind speed u, potential temperature \theta.
* \mathbf{y} (Fast variables): Turbulent kinetic energy e, Reynolds stresses \overline{u'w'}.
* \epsilon: Ratio of turbulent to synoptic timescales.

2. Architectural Overview

The architecture is partitioned into seven functional domains, ensuring a clean separation between observation, geometry, and numerical discretization.

Table 1: Framework Subsystems

Layer	Primary Submodule	Scientific Role
Observation	DataIngestion.jl	Ingests CSV/NetCDF data; projects profiles onto Gegenbauer polynomials.
Geometry	FoldTracking.jl / CanardDetection.jl	Tracks slow manifolds (\mathcal{S}_0); detects fold points and normal hyperbolicity.
Discovery	WSINDyEngine.jl / ModelSelection.jl	Constructs weak-form libraries; evaluates AIC/BIC Pareto trade-off curves.
Closures	MOSTClosure.jl / WSINDyClosure.jl	Houses C^\infty-differentiable eddy diffusivity profiles (K_m, K_h).
Discretization	StretchedGrid.jl / SpectralBLGalerkin.jl	Offers vertical finite differences and zero-allocation spectral expansions.
Calibration	HierarchicalTuring.jl	Performs Bayesian MCMC sampling and multi-site hierarchical pooling.
Reporting	LaTeXExporter.jl	Converts symbolic ASTs into publication-ready LaTeX equations.

3. Mathematical Foundations

3.1 GSPT and Turbulence Collapse

In the singular limit (\epsilon \to 0), transients decay onto a critical manifold \mathcal{S}_0. By Fenichel’s Theorem, this manifold persists as long as it is normally hyperbolic, meaning the eigenvalues of the fast Jacobian \mathbf{D}_{\mathbf{y}}\mathbf{g} have non-zero real parts.

Turbulence collapse occurs at the fold locus, where normal hyperbolicity is lost: \det\left(\mathbf{D}_{\mathbf{y}}\mathbf{g}(\mathbf{z})\right) = 0

ASM.jl applies adjugate desingularization to integrate trajectories across these singular lines, allowing for the detection of canard trajectories—solutions that follow unstable manifold sheets without immediate turbulent breakdown.

3.2 Weak-Form WSINDy Discovery

The discovery engine avoids noise by projecting PDEs against C^\infty test functions \phi_{i,j}(z, t): \iint_{\Omega} u(z, t) \, \frac{\partial \phi_{i,j}}{\partial t} \, dz \, dt + \iint_{\Omega} K_m(z) \frac{\partial u}{\partial z} \, \frac{\partial \phi_{i,j}}{\partial z} \, dz \, dt = 0

Physical Hard Constraints

To ensure realism, candidate coefficients are solved via constrained quadratic programming, enforcing:

1. Positivity: K_m(z) \ge 0 and K_h(z) \ge 0.
2. Energy Monotonicity: Non-negative column integration of dissipation terms.
3. Neutral Limit Recovery: Convergence to MOST behavior as Ri \to 0.

3.3 C^\infty Regularity

Numerical continuation engines (like BifurcationKit.jl) require continuous Jacobians. ASM.jl replaces piecewise "kinks" with smooth algebraic operators:

* Smooth Max: \operatorname{smooth\_max}(x, y; \epsilon) = \frac{x + y + \sqrt{(x - y)^2 + \epsilon^2}}{2}
* Regime Blending: Hyperbolic tangent functions blend stable and unstable regimes smoothly: w(\zeta) = \frac{1}{2} [ 1 + \tanh(\zeta/\epsilon) ].

4. Discretization and Numerical Efficiency

ASM.jl provides two primary spatial backends:

1. MethodOfLinesFD: Uses stretched 1D vertical finite differences. A stretching parameter \alpha concentrates grid points near the surface where shear and thermal gradients are steepest.
2. SpectralBLGalerkin: Expands variables in orthogonal Gegenbauer polynomials C_n^{(\lambda)}(z).
  * Zero-Allocation: Evaluation of non-linear products is performed in modal space via precomputed 3-tensor contraction matrices C_{ijk}^{(\lambda)}.
  * This eliminates dynamic memory allocation (0 Bytes) in stiff ODE loops (e.g., Rodas5P, RadauIIA5).

5. Calibration and Uncertainty Quantification

5.1 Model Selection

Candidate models are evaluated along a Pareto front comparing complexity (number of active terms k) against error (RSS or R^2). Selection is guided by:

* AIC: 2k + N \ln(\text{RSS}/N)
* BIC: k \ln(N) + N \ln(\text{RSS}/N)

5.2 Hierarchical Bayesian Inference

HierarchicalTuring.jl pools statistical strength across distinct field sites.

* Global Hyperpopulation: Represents universal constants (e.g., von Kármán constant \kappa).
* Site-Specific Parameters: Adapt to local roughness (z_0) or soil moisture while remaining constrained by the global mean.

5.3 Uncertainty Propagation

The system forwards parameter posterior chains into spatial distributions of eddy diffusivity (K_m), providing 95% credibility ribbons:

* Median: Q_{0.50}
* Bounds: Q_{0.025} and Q_{0.975}

6. Observational Benchmarks

The framework has been validated against four major field campaigns.

Table 2: Field Campaign Diagnostics

Campaign	Observations	Mean Stability (\bar{\zeta})	Observed Dynamics
FLOSS	70,796	0.3842	Extreme manifold contraction over ice/snow; rapid mode damping.
CASES-99	6,538	0.2114	Strong non-equilibrium transients during sunset collapse.
BLLAST	5,600	0.1985	Active fast-slow energy exchange during evening transitions.
SHEBA	2,273	0.3450	Strong manifold alignment in stable Arctic layers.

7. Key Software Capabilities

* Symbolic Interoperability: Discovered models export directly to ModelingToolkit.jl, enabling automated generation of analytical Jacobians.
* Automated LaTeX Pipeline: The LaTeXExporter.jl module converts discovered symbolic terms into publication-ready differential equations (e.g., \frac{\partial u}{\partial t} = \nu \frac{\partial^2 u}{\partial z^2} + \xi_1 z \frac{\partial u}{\partial z}).
* Comprehensive Diagnostics: The Diagnostics module tracks statistical metrics (RMSE, R^2), physical energy budgets (TKE), and geometric manifold metrics (fold distance, transversality).

"Capturing these [canard] trajectories allows ASM.jl to prevent spurious numerical oscillations during transient nocturnal transitions." — Source Analysis on Canard Detection.
