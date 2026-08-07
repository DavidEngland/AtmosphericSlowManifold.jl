`ErrorMetrics` provides high-performance, allocation-minimal evaluation metrics for comparing predicted atmospheric profiles and learned closures ($y_{\text{pred}}$) against observational data or reference trajectories ($y_{\text{true}}$).

---

### Mathematical Metrics Summary

| Function | Mathematical Formulation | Division by Zero Protection |
| --- | --- | --- |
| **`rmse`** | $\text{RMSE} = \sqrt{\frac{1}{N}\sum_{i=1}^N (y_{\text{pred},i} - y_{\text{true},i})^2}$ | Returns `0.0` for empty arrays. |
| **`mae`** | $\text{MAE} = \frac{1}{N}\sum_{i=1}^N \Vert{}y_{\text{pred},i} - y_{\text{true},i}\Vert{}$ | Returns `0.0` for empty arrays. |
| **`bias`** | $\text{Bias} = \frac{1}{N}\sum_{i=1}^N (y_{\text{pred},i} - y_{\text{true},i})$ | Returns `0.0` for empty arrays. |
| **`r2`** | $R^2 = 1 - \frac{\sum (y_{\text{true},i} - y_{\text{pred},i})^2}{\sum (y_{\text{true},i} - \bar{y}_{\text{true}})^2}$ | Returns `1.0` if $\text{Var}(y_{\text{true}}) = 0$. |
| **`nrmse`** | $\text{NRMSE} = \frac{\text{RMSE}}{S}, \quad S \in \{\sigma(y_{\text{true}}), \Delta y_{\text{true}}, \Vert{}\bar{y}_{\text{true}}\Vert{}\}$ | Returns raw $\text{RMSE}$ if scale $S = 0$. |
| **`skill_score`** | $SS = 1 - \frac{\text{MSE}_{\text{pred}}}{\text{MSE}_{\text{ref}}} = 1 - \frac{\sum (y_{\text{pred},i} - y_{\text{true},i})^2}{\sum (y_{\text{ref},i} - y_{\text{true},i})^2}$ | Returns `0.0` if $\text{MSE}_{\text{ref}} = 0$. |
| **`correlation`** | $r = \frac{\sum (y_{\text{pred},i} - \bar{y}_{\text{pred}})(y_{\text{true},i} - \bar{y}_{\text{true}})}{\sqrt{\sum (y_{\text{pred},i} - \bar{y}_{\text{pred}})^2 \sum (y_{\text{true},i} - \bar{y}_{\text{true}})^2}}$ | Returns `0.0` if standard deviation is zero. |
| **`normalized_bias`** | $\text{NBias} = \frac{\text{Bias}}{\Vert{}\bar{y}_{\text{true}}\Vert{}}$ | Fallbacks to unnormalized `bias` if $\bar{y}_{\text{true}} = 0$. |
| **`relative_l2_error`** | $\text{Rel } L_2 = \frac{\Vert{}y_{\text{pred}} - y_{\text{true}}\Vert{}_2}{\Vert{}y_{\text{true}}\Vert{}_2} = \sqrt{\frac{\sum (y_{\text{pred},i} - y_{\text{true},i})^2}{\sum y_{\text{true},i}^2}}$ | Returns absolute $L_2$ error if $\Vert{}y_{\text{true}}\Vert{}_2 = 0$. |
| **`closure_residual`** | Pointwise residual vector: $r_i = \Vert{} \phi_{\text{obs},i} - \phi_{\text{model},i} \Vert{}$ | Allocates destination output array matching input dimension. |

---

### Key Implementation Features

1. **Inbounds Vectorization:** Loop operations use `@inbounds for i in eachindex(...)` to skip bounds checking during trajectory scans.
2. **Type-Safe `Float64` Accumulation:** Explicit casting (`Float64(...)`) inside accumulator loops prevents integer overflow or single-precision float accuracy degradation during accumulation over large grids.
3. **Array Dimension Assertion:** `_assert_same_length` enforces matching dimensions across $y_{\text{pred}}$, $y_{\text{true}}$, and $y_{\text{ref}}$ prior to entering computation loops, throwing clear `ArgumentError` exceptions on dimension mismatches.