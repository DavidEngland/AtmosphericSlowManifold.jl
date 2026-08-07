`src/Manifold/GSPTDiagnostics.jl` provides diagnostic routines to quantify numerical proximity to manifold fold singularities ($\det(D_{\mathbf{y}}\mathbf{g}) = 0$) and evaluate normal transversality conditions using symbolic differentiation via `Symbolics.jl`.

### Structure & Function Breakdown

* **`FoldConstraint`**: Encapsulates the symbolic algebraic condition $\mathcal{F}(\mathbf{z}) = 0$ that defines the fold manifold where normal hyperbolicity degrades.
* **`fold_residual(fold, substitutions)`**: Substitutes real-valued state variables into `fold.expr` and extracts the numeric scalar residual using `Symbolics.value`. A residual near zero signals proximity to a turbulence collapse threshold or fold bifurcation point.
* **`fold_transversality(fold, variable)`**: Computes the exact symbolic derivative $\frac{\partial \mathcal{F}}{\partial z_i}$ of the fold expression with respect to a target coordinate $z_i$.

---

### Execution & Diagnostic Role

```
       Physical Data / State
                 │
                 ▼
     [ substitutions Dict ]
                 │
                 ├───> fold_residual() ──────> Proximity Metric: F(z) ≈ 0
                 │
                 └───> fold_transversality() ─> Gradient / Transversality: ∂F/∂z_i

```

1. **Observational Validation:** Used during field data evaluation (e.g., FLOSS, CASES-99) to score individual observational time steps based on how closely they align with the singular locus.
2. **Bifurcation Tracking:** Interoperates with `Geometry/` and `BifurcationKit.jl` by providing exact analytical gradients required for continuation algorithms when tracing fold curves and Hopf points.