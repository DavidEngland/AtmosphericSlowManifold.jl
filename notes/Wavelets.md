Integrating wavelet-based wave-turbulence decomposition directly addresses Dick's core critique—proving that high-$Ri$ states reflect active manifold dynamics rather than passive wave variance. Continuous Morlet wavelets give you the exact time-frequency resolution needed to isolate intermittent turbulence bursts from non-transporting gravity wave packets in raw 20 Hz sonic data.

**Preprocessing & Wavelet Pipeline Architecture**

* **Continuous Morlet Wavelet Transform (CWT):** Apply CWT to 20 Hz $u, v, w, T$ time series. Morlet wavelets preserve phase information, making them ideal for detecting localized, non-stationary wave packets versus coherent turbulent eddies.
* **Cospectral Gap Identification:** Identify the dynamic scale-dependent minimum in the real part of the cross-wavelet spectrum (cospectra of $w'\theta'$ and $u'w'$). Scales larger than the gap represent non-transporting wave motion or wind-meandering; smaller scales represent active turbulent flux.
* **Variance Partitioning:** Separate total variance into turbulent kinetic energy ($e_{turb}$) and wave kinetic energy ($e_{wave}$). Recompute gradient Richardson numbers using wave-filtered turbulent scales to produce a stabilized $Ri_{clean}$.

**Julia Workflow (`AtmosphericSlowManifold.jl`)**

| Analysis Stage | Tool / Method | Operational Goal |
| --- | --- | --- |
| **Data Ingestion** | `DSP.jl` / Custom I/O | Load 20 Hz sonic anemometer series ($u, v, w, T$) from CASES-99 / SHEBA |
| **Multiscale CWT** | `ContinuousWavelets.jl` | Run Morlet transform to construct scale-dependent cross-power spectra |
| **Gap Extraction** | Cospectral Gap Algorithm | Partition scalar and momentum fluxes where cospectral density drops to zero |
| **Fold Sharpening** | Manifold Diagnostic Engine | Re-plot $Ri_{clean}$ against non-dimensional gradients to isolate kinematic fold "knees" |

**Key Diagnostic Outcomes**

* **$Ri$ Noise Stabilization:** Removing $e_{wave}$ eliminates spurious high-$Ri$ spikes caused by non-mixing orbital wave velocities.
* **Fold Sharpening:** Stripping wave variance removes the blurring around transition "knees," clarifying whether the fold structure represents a true physical bifurcation on the slow manifold.
* **Automated Quality Control:** Integrates cleanly with SHEBA's FFT quality flags (`fl1`–`fl5`), providing a reproducible, automated QC pipeline before feeding profiles into your manifold solver.

Would you like to start by sketching the CWT cospectral gap extraction module in Julia, or focus on setting up the pipeline to ingest the CASES-99 20 Hz dataset first?