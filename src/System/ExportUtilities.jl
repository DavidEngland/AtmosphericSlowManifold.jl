module ExportUtilities

using CSV
using DataFrames
using JSON3
using NCDatasets

using ..AtmosphericSlowManifold: ObservationTable
using ..AtmosphericSlowManifold: ModalBudgetDiagnostic
using ..AtmosphericSlowManifold: DiscoveredModel
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

function _json_ready(x)
    if x === nothing || x isa AbstractString || x isa Bool
        return x
    elseif x isa Number
        xf = Float64(x)
        return isfinite(xf) ? x : nothing
    elseif x isa Symbol
        return String(x)
    elseif x isa AbstractVector
        return [_json_ready(v) for v in x]
    elseif x isa Dict
        out = Dict{String, Any}()
        for (k, v) in x
            out[String(k)] = _json_ready(v)
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
