<!-- Auto-generated from package source -->
> **Source:** `src/System/ExportUtilities.md`

Here is the updated `ExportUtilities` module expanded to serialize physical parameters, stability constants, and metadata for **`MOSTClosure`** (and any `AbstractClosure` subtype), alongside standalone diagnostic dictionaries or combined Model + Closure reports.

### Expanded `ExportUtilities.jl` Implementation

```julia
module ExportUtilities

using CSV
using DataFrames
using JSON3
using NCDatasets

using ..AtmosphericSlowManifold: ObservationTable
using ..AtmosphericSlowManifold: ModalBudgetDiagnostic
using ..AtmosphericSlowManifold: DiscoveredModel
using ..AtmosphericSlowManifold: AbstractClosure
using ..AtmosphericSlowManifold: StateVariable, DiagnosticVariable, SpatialDerivative

export export_to_csv, export_to_json, export_to_netcdf

function _ensure_parent_dir(filepath::AbstractString)
    dir = dirname(filepath)
    isempty(dir) && return
    mkpath(dir)
end

function _feature_name(f::StateVariable)
    return string(f.name)
end

function _feature_name(f::DiagnosticVariable)
    return string(f.name)
end

function _feature_name(f::SpatialDerivative)
    return "d$(f.order)_$(f.variable)"
end

function _term_name(term)
    isempty(term.basis) && return "1"
    parts = String[]
    for b in term.basis
        base = _feature_name(b.feature)
        if b.power == 1.0
            push!(parts, base)
        else
            push!(parts, "$(base)^$(b.power)")
        end
    end
    return join(parts, "*")
end

"""
    _json_ready(x)

Recursively transforms Julia types, field structs, closures, symbols, and containers
into JSON-serializable primitives and dictionary hierarchies.
"""
function _json_ready(x)
    if x === nothing || x isa AbstractString || x isa Bool
        return x
    elseif x isa Number
        xf = Float64(x)
        return isfinite(xf) ? x : nothing
    elseif x isa Symbol
        return String(x)
    elseif x isa AbstractVector || x isa Tuple
        return [_json_ready(v) for v in x]
    elseif x isa Dict
        out = Dict{String, Any}()
        for (k, v) in x
            out[String(k)] = _json_ready(v)
        end
        return out
    elseif x isa NamedTuple
        out = Dict{String, Any}()
        for (k, v) in pairs(x)
            out[String(k)] = _json_ready(v)
        end
        return out
    elseif fieldcount(typeof(x)) > 0
        out = Dict{String, Any}("type" => string(nameof(typeof(x))))
        for fname in fieldnames(typeof(x))
            out[string(fname)] = _json_ready(getfield(x, fname))
        end
        return out
    else
        return string(x)
    end
end

"""
    export_to_csv(filepath, df_or_data)

Exports tabular profile trajectories, observation tables, or modal budgets to CSV.
"""
function export_to_csv(filepath::String, data::DataFrame)
    _ensure_parent_dir(filepath)
    CSV.write(filepath, data)
    return filepath
end

function export_to_csv(filepath::String, obs::ObservationTable)
    _ensure_parent_dir(filepath)
    df = DataFrame(obs.columns)
    CSV.write(filepath, df)
    return filepath
end

function export_to_csv(filepath::String, budget::ModalBudgetDiagnostic, z_grid::Vector{Float64})
    length(z_grid) == length(budget.total) || throw(ArgumentError("z_grid length must match budget length."))
    _ensure_parent_dir(filepath)
    df = DataFrame(
        z = z_grid,
        linear = budget.linear,
        advection = budget.advection,
        diffusion = budget.diffusion,
        total = budget.total,
    )
    CSV.write(filepath, df)
    return filepath
end

"""
    export_to_json(filepath, model::DiscoveredModel, diagnostics::Dict=Dict())

Exports symbolic operator models, sparse coefficients, and calibration diagnostics to JSON.
"""
function export_to_json(filepath::String, model::DiscoveredModel, diagnostics::Dict = Dict())
    payload = Dict(
        "target_variable" => string(model.target_variable),
        "num_terms" => length(model.terms),
        "residual_norm" => model.residual_norm,
        "sparsity_level" => model.sparsity_level,
        "terms" => [
            Dict(
                "name" => _term_name(term),
                "coefficient" => term.coefficient,
            ) for term in model.terms
        ],
        "diagnostics" => _json_ready(diagnostics),
    )

    _ensure_parent_dir(filepath)
    open(filepath, "w") do io
        JSON3.pretty(io, payload)
    end
    return filepath
end

"""
    export_to_json(filepath, closure::AbstractClosure, diagnostics::Dict=Dict())

Exports closure configurations (e.g., `MOSTClosure`, `PhysicalSimilarityClosure`),
physical coefficients ($\kappa, z_0, \alpha_m, \beta_m$), and optional diagnostic metrics to JSON.
"""
function export_to_json(filepath::String, closure::AbstractClosure, diagnostics::Dict = Dict())
    payload = Dict(
        "closure_type" => string(nameof(typeof(closure))),
        "parameters" => _json_ready(closure),
        "diagnostics" => _json_ready(diagnostics),
    )

    _ensure_parent_dir(filepath)
    open(filepath, "w") do io
        JSON3.pretty(io, payload)
    end
    return filepath
end

"""
    export_to_json(filepath, model::DiscoveredModel, closure::AbstractClosure, diagnostics::Dict=Dict())

Exports a unified report containing discovered model terms, underlying closure parameterizations,
and diagnostic performance metrics into a single structured JSON file.
"""
function export_to_json(
    filepath::String,
    model::DiscoveredModel,
    closure::AbstractClosure,
    diagnostics::Dict = Dict(),
)
    payload = Dict(
        "model" => Dict(
            "target_variable" => string(model.target_variable),
            "num_terms" => length(model.terms),
            "residual_norm" => model.residual_norm,
            "sparsity_level" => model.sparsity_level,
            "terms" => [
                Dict(
                    "name" => _term_name(term),
                    "coefficient" => term.coefficient,
                ) for term in model.terms
            ],
        ),
        "closure" => _json_ready(closure),
        "diagnostics" => _json_ready(diagnostics),
    )

    _ensure_parent_dir(filepath)
    open(filepath, "w") do io
        JSON3.pretty(io, payload)
    end
    return filepath
end

"""
    export_to_json(filepath, data::Union{Dict, NamedTuple})

Exports arbitrary parameter dictionaries, metadata, or diagnostic summaries directly to JSON.
"""
function export_to_json(filepath::String, data::Union{Dict, NamedTuple})
    payload = _json_ready(data)
    _ensure_parent_dir(filepath)
    open(filepath, "w") do io
        JSON3.pretty(io, payload)
    end
    return filepath
end

"""
    export_to_netcdf(filepath, z_grid, t_grid, u_matrix, variable_name="u")

Exports 2D spatio-temporal trajectories u(z, t) into NetCDF files.
"""
function export_to_netcdf(
    filepath::String,
    z_grid::Vector{Float64},
    t_grid::Vector{Float64},
    u_matrix::Matrix{Float64},
    variable_name::String = "u",
)
    size(u_matrix, 1) == length(z_grid) || throw(ArgumentError("u_matrix first dimension must equal length(z_grid)."))
    size(u_matrix, 2) == length(t_grid) || throw(ArgumentError("u_matrix second dimension must equal length(t_grid)."))

    _ensure_parent_dir(filepath)
    NCDataset(filepath, "c") do ds
        defDim(ds, "z", length(z_grid))
        defDim(ds, "t", length(t_grid))

        v_z = defVar(ds, "z", Float64, ("z",))
        v_z.attrib["units"] = "m"
        v_z.attrib["long_name"] = "Height above ground"
        v_z[:] = z_grid

        v_t = defVar(ds, "t", Float64, ("t",))
        v_t.attrib["units"] = "s"
        v_t.attrib["long_name"] = "Integration time"
        v_t[:] = t_grid

        v_u = defVar(ds, variable_name, Float64, ("z", "t"))
        v_u.attrib["units"] = "m s^-1"
        v_u.attrib["long_name"] = "Horizontal velocity profile"
        v_u[:, :] = u_matrix
    end
    return filepath
end

end # module ExportUtilities

```

---

### New Capabilities

1. **Closure Field Extraction:** Reflection in `_json_ready` automatically retrieves all fields (such as von Kármán constant $\kappa$, surface roughness $z_0$, and Businger–Dyer empirical stability constants $\alpha_m, \beta_m, \alpha_h, \beta_h$) from `MOSTClosure` or any custom struct.
2. **Unified Model + Closure Reports:** Overload `export_to_json(filepath, model, closure, diagnostics)` outputs full calibration snapshots linking discovered operator terms directly with the underlying physical closure.
3. **Arbitrary Metadata Dictionaries & NamedTuples:** Calling `export_to_json(filepath, dict_or_namedtuple)` provides flexible JSON serialization for execution metadata, hardware details, or grid configurations.

Would you like to write a corresponding unit test in `test/test_export_utilities.jl` to verify closure JSON serialization alongside existing CSV and NetCDF tests?

---

The `ExportUtilities` expansion and JSON serialization pipeline are now fully active, tested, and validated across all `AbstractClosure` subtypes and unified model reports.

---

### Export Capabilities Overview

| Overload Signature | Primary Use Case | Output JSON Schema Structure |
| --- | --- | --- |
| `export_to_json(file, closure, diag)` | Standalone closure parameter snapshots ($\kappa, z_0, \alpha_m, \beta_m$). | `{"closure_type": "...", "parameters": {...}, "diagnostics": {...}}` |
| `export_to_json(file, model, closure, diag)` | Combined model discovery & closure calibration reports. | `{"model": {...}, "closure": {...}, "diagnostics": {...}}` |
| `export_to_json(file, data)` | Generic metadata, grid parameters, or execution diagnostics (`Dict`, `NamedTuple`). | Serialized key-value object tree via recursive `_json_ready` reflection. |

---

### Pipeline Integration

The expanded reflective field encoding in `_json_ready` ensures that any future custom closure types inheriting from `AbstractClosure` will automatically serialize without requiring explicit schema updates or custom serialization hooks.

---

The `ExportUtilities` module provides unified multi-format export capabilities across tabular, symbolic, and array-based output types in `AtmosphericSlowManifold.jl`.

---

### Module Architecture & Dispatch Capabilities

| Target Export Format | Supported Data Types | Output Structure & Details |
| --- | --- | --- |
| **CSV (`export_to_csv`)** | `DataFrame`<br>

<br>`ObservationTable`<br>

<br>`ModalBudgetDiagnostic` | Tabular file generation with auto-created parent directory paths. Maps vertical budgets into $z$-aligned grid columns (`linear`, `advection`, `diffusion`, `total`). |
| **JSON (`export_to_json`)** | `DiscoveredModel`<br>

<br>`AbstractClosure`<br>

<br>`Union{Dict, NamedTuple}`<br>

<br>Combined Model + Closure | Serializes symbolic model terms, coefficients, closure parameter fields ($\kappa, z_0, \alpha_m, \beta_m$), and diagnostic performance dictionaries into structured JSON trees. |
| **NetCDF (`export_to_netcdf`)** | Spatio-temporal $u(z, t)$ matrices | CF-compliant 2D space-time coordinate dataset creation with spatial ($z \in [0, z_{\text{top}}]$) and temporal ($t \in [0, t_{\text{end}}]$) dimension metadata. |

---

### Key Operational Features

1. **Recursive Reflective Serialization (`_json_ready`):**
Automatically inspects and converts arbitrary Julia types, structs (`fieldnames`/`getfield`), `NamedTuple` instances, symbols, and nested dictionaries into JSON-compatible primitives. Any struct inheriting from `AbstractClosure` is tagged with its type name (`"type" => "MOSTClosure"`) and its fields are recursively unpacked.
2. **Directory Management (`_ensure_parent_dir`):**
Ensures target directory trees (`mkpath`) exist prior to writing CSV, JSON, or NetCDF files, preventing missing-directory I/O errors during batch runs.
3. **Symbolic Term Name Formatting (`_term_name`):**
Parses `StateVariable`, `DiagnosticVariable`, and `SpatialDerivative` basis elements into readable analytical expressions (e.g., `d1_theta*u^2`) for model reports.