`scripts/generate_report_supplements.jl` is a Julia utility script that generates LaTeX and Markdown supplement tables for repository reproducibility, numerical verification, and cross-campaign benchmark results.

---

### Core Functionality & Output Artifacts

#### 1. Reproducibility Metadata (`reproducibility.tex` & `reproducibility.md`)

* Extracts repository state via `git` (`remote url`, `short commit`).
* Inspects system environment (`Julia version`, `threads`, `kernel`, `CPU architecture`).
* Computes `Manifest.toml` SHA-256 hash digest (first 16 characters).
* Documents uncertainty quantification settings ($N = 1,000$ draws, `seed = 42`).

#### 2. Numerical Verification Checks (`numerical_verification.tex` & `numerical_verification.md`)

* **Solver Integrations:** Verifies all physical and baseline ODE/PDE solvers return `"Success"`.
* **Metric Identical Assertions:** Runs `candidate_metrics_valid()` to verify:

$$\text{RSS} \approx \text{residual\_norm}^2$$


$$\text{AIC} \approx 2k + N \ln\left(\frac{\text{RSS}}{N}\right)$$


$$\text{BIC} \approx k \ln(N) + N \ln\left(\frac{\text{RSS}}{N}\right)$$


* **Pareto Non-Emptiness:** Ensures non-zero candidate selection across both $RSS$ and $R^2$ fronts.

#### 3. Cross-Campaign Benchmark Summary (`cross_campaign_results.tex` & `cross_campaign_results.md`)

* Parses `benchmark_summary.json` across observation sources (e.g., CASES-99, SHEBA).
* Formats best-model candidate selections along with $R^2$, residual norm, $AIC$, and $BIC$.