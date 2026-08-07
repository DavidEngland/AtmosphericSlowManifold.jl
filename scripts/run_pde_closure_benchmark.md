The `scripts/run_pde_closure_benchmark.jl` script executes an end-to-end PDE simulation, data-driven closure selection, and LaTeX/figure artifact generation pipeline across field campaign datasets.

---

## 1. Architectural & Data Flow Overview

The script connects preprocessed field observations (`reports/generated/campaign_exports/`) with the core `AtmosphericSlowManifold.jl` solving and model-selection modules:

```
┌──────────────────────────────────────────┐
┌ Campaign Exports Data                   │
│ • JSON: campaign model & diagnostics     │
│ • CSV: raw observation profiles          │
└────────────────────┬─────────────────────┘
                     │
                     ▼
┌──────────────────────────────────────────┐
│ Closure System Construction              │
│ • PhysicalSimilarityClosure (Discovered) │
│ • PhysicalSimilarityClosure (Baseline)   │
└────────────────────┬─────────────────────┘
                     │
                     ▼
┌──────────────────────────────────────────┐
│ Galerkin Discretization & ODE Solve      │
│ • SpectralBLGalerkin (N = 12 modes)      │
│ • Fallback Chain: Rodas5P → RadauIIA5    │
└────────────────────┬─────────────────────┘
                     │
                     ▼
┌──────────────────────────────────────────┐
│ Pareto Model Selection & Export          │
│ • Evaluate RSS, R², AIC, BIC             │
│ • Compute non-dominated Pareto sets      │
│ • Export LaTeX, PNG, CSV, and JSON3      │
└──────────────────────────────────────────┘

```

---

## 2. Key Subsystem Implementation Breakdown

### 2.1 Stiff Solver Fallback Chain (`solve_with_fallback`)

To handle stiff ODE systems resulting from high-order Gegenbauer spectral Galerkin discretizations ($N = 12$), the solver iterates through a sequence of stiffly stable implicit integrators:

1. **`Rodas5P()`**: 5th-order A-stable stiffly-accurate Rosenbrock-Wanner scheme with adaptive error control.
2. **`RadauIIA5()`**: 5th-order fully implicit Runge-Kutta method for extremely stiff boundary-layer transitions.
3. **`FBDF()`**: Fixed-leading-coefficient Backward Differentiation Formula for high-dimensional stiff systems.

### 2.2 Pareto Selection & Multi-Objective Optimization

Candidate models (baseline vs. physically constrained discovered closures) are evaluated across both residual sum of squares ($\text{RSS}$) and variance explained ($R^2$):

$$\text{AIC} = 2k + N \ln\left(\frac{\text{RSS}}{N}\right), \qquad \text{BIC} = k \ln(N) + N \ln\left(\frac{\text{RSS}}{N}\right)$$

Non-dominated solutions are extracted using `compute_pareto_front` to isolate optimal complexity-versus-accuracy balances.

---

## 3. Produced Output Artifact Matrix

Executing `julia --project=. scripts/run_pde_closure_benchmark.jl` populates `reports/generated/pde_benchmark/` with the following assets:

| Directory | Artifact File | Type | Scientific Role / Description |
| --- | --- | --- | --- |
| `tables/` | `pde_benchmark_summary.csv` | CSV | Tabular metric comparison ($\text{RMSE}$, $R^2$, $\text{AIC}$, $\text{BIC}$, $K_m^{\max}$) across campaign runs. |
| `tables/` | `cross_campaign_best_models.tex` | LaTeX | Production-ready cross-campaign summary table generated via `latex_site_summary_table`. |
| `tables/` | `<site>_best_equation.tex` | LaTeX | Display equation block ($\\[\partial_t u = \dots \\]$) for the top Pareto closure. |
| `tables/` | `<site>_best_terms.tex` | LaTeX | Individual term coefficient table with standard errors and residual norms. |
| `figures/` | `<site>_pareto_rss.png` | Plot | Model complexity $k$ vs. $\text{RSS}$ trade-off curve with Pareto front highlighting. |
| `figures/` | `<site>_pareto_r2.png` | Plot | Model complexity $k$ vs. $R^2$ accuracy curve with Pareto front highlighting. |
| `figures/` | `pde_profile_comparison.png` | Plot | Multi-panel modal profile trajectory comparisons (Observed Initial vs. Baseline vs. Discovered). |
| `./` | `benchmark_summary.json` | JSON3 | Full machine-readable manifest containing exact candidate counts, LaTeX strings, and performance metrics. |

---

## 4. Verification & Execution Target

To run the complete pipeline and re-generate all benchmark deliverables in the project workspace:

```bash
make pde-benchmark

```