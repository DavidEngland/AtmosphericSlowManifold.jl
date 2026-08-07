`ObservationTable` and its accompanying ingestion pipeline provide flexible data loading, unit verification, and header normalization for tower, radiosonde, and LES profile observations across CSV and NetCDF formats.

---

### Ingestion Pipeline Architecture

```
        Raw CSV / NetCDF Files (.csv, .nc)
                       │
                       ▼
       Header & Column Alias Mapping
   (HEADER_ALIASES, SURFACE_FLUX_ALIASES)
                       │
                       ▼
       Type Conversion & Safe Parsing
          (_to_f64, _normalize_table)
                       │
                       ▼
        Derived Physics Ingestion
      (_derive_obukhov_length: L)
                       │
                       ▼
     Strict / Expected Unit Assertions
   (EXPECTED_UNITS, _assert_units)
                       │
                       ▼
               ObservationTable
   (Dict{Symbol, Vector{Float64}}, Units)

```

---

### Core Data Structures & Constants

| Component | Target / Type | Function / Specification |
| --- | --- | --- |
| `ObservationTable` | Struct | In-memory container mapping column symbols to `Vector{Float64}` data arrays alongside unit metadata (`Dict{Symbol, String}`). |
| `REQUIRED_TOWER_COLUMNS` | Constant Vector | Ensures required fields exist: `:z`, `:u`, `:v`, `:theta`, `:q`, `:u_star`. |
| `HEADER_ALIASES` | Constant Dict | Maps canonical column symbols to candidate header aliases (e.g., `:theta` $\to$ `:theta_k`, `:potential_temperature_k`, `:theta`). |
| `EXPECTED_UNITS` | Constant Dict | Canonical units: $z$ (`m`), $u, v, u_*$ (`m s^-1`), $\theta$ (`K`), $q$ (`kg kg^-1`). |

---

### Physical Formulations

#### Monin–Obukhov Length Scale ($L$)

When sensible heat flux ($H_s$) and friction velocity ($u_*$) are present, `_derive_obukhov_length` computes the Obukhov length profile:

$$L = -\frac{\theta_0 \cdot u_*^3 \cdot \rho c_p}{\kappa \cdot H_s}$$

* $\theta_0$: Reference potential temperature ($\text{K}$).
* $u_*:$ Friction velocity ($\text{m/s}$).
* $\rho c_p$: Volumetric heat capacity of air ($\approx 1200.0 \text{ J m}^{-3} \text{ K}^{-1}$).
* $\kappa$: von Kármán constant ($0.40$).
* $H_s$: Sensible heat flux ($\text{W m}^{-2}$).
* Safeguards against division-by-zero singularities when $\vert{}H_s\vert{} \le 10^{-12} \text{ W m}^{-2}$ by returning `NaN`.

---

### Key Capabilities & Utility API

* **Unified Entry Point (`read_observation_data`):** Single function call handling both `.csv` and `.nc` file extensions with column selection fallbacks and default values for missing variables (e.g., $q = 0.0$, $u_* = 0.3 \text{ m/s}$).
* **Strict & Flexible Unit Verification:** Supports strict header unit suffixes (`_m`, `_k`, `_ms`, `_kgkg`) or automatic assignment to expected SI standards.
* **Sibling Directory Resolver (`resolve_sibling_data_dir`):** Locates dataset paths relative to sibling repositories (such as `SpectralBL-Analytics`) within workspace environments.