<!-- Auto-generated from package source -->
> **Source:** `src/Discretization/Backends/MethodOfLinesFD.md`

`src/Discretization/Backends/MethodOfLinesFD.jl` implements a concrete spatial discretization backend for `PDESystem` handling via `MethodOfLines.jl`. It maps continuous boundary layer equations onto stretched 1D vertical grids and compiles them into stiff differential equation problems for `DifferentialEquations.jl`.

---

### Key Operational Characteristics

1. **Grid Stretching Parameterization ($\alpha$):**
The `generate_stretched_grid` helper uses parameters $N$ (node count), $H$ (domain height $z_{\text{top}}$), and $\alpha$ (stretching exponent) to concentrate vertical grid resolution near the surface boundary layer ($z \to 0$) where turbulent shear gradients are steepest:

$$z_i = H \left( \frac{i}{N} \right)^\alpha, \qquad i \in \{0, 1, \dots, N\}$$


2. **`MethodOfLines.jl` Compilation (`MOLFiniteDifference`):**
`MOLFiniteDifference` maps spatial derivatives ($\partial_z, \partial_{zz}$) along independent variable `zvar` using finite differences of specified order (`disc.order`), transforming spatial PDE fields into a system of coupled ODEs/DAEs in time `tvar`.
3. **Stiff Time Integration Default (`TRBDF2`):**
Defaulting `solver = TRBDF2()` provides implicit, $L$-stable time integration optimized for stiff parabolic diffusion systems and non-linear boundary flux interactions.

---

### Backend Integration Summary

| Field / Parameter | Default | Physical & Numerical Purpose |
| --- | --- | --- |
| `N::Int` | `100` | Number of vertical grid intervals across the domain. |
| `H::Float64` | `3000.0` | Boundary layer domain height $z_{\text{top}}$ (in meters). |
| `alpha::Float64` | `3.5` | Surface mesh clustering exponent ($\alpha > 1$ refines grid near $z = 0$). |
| `order::Int` | `2` | Spatial finite difference approximation order. |
| `solver` | `TRBDF2()` | Stiff implicit ODE integrator with automatic step-size control. |