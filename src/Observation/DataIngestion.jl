# AtmosphericSlowManifold.jl: Data ingestion utilities for observational tower data
# src/Observation/DataIngestion.jl
using DelimitedFiles
using NCDatasets
using CSV
using DataFrames

struct ObservationTable
    columns::Dict{Symbol,Vector{Float64}}
    units::Dict{Symbol,String}
end

const REQUIRED_TOWER_COLUMNS = [:z, :u, :v, :theta, :q, :u_star]

const HEADER_ALIASES = Dict(
    :z => (:z_m, :height_m, :z, :height, :z_5m),
    :u => (:u_ms, :u_m_s, :u, :u_wind, :u_component),
    :v => (:v_ms, :v_m_s, :v, :v_wind, :v_component),
    :theta => (:theta_k, :potential_temperature_k, :theta, :potential_temp, :potential_temperature),
    :q => (:q_kgkg, :specific_humidity_kgkg, :q, :specific_humidity),
    :u_star => (
        :u_star_ms, :ustar_ms, :friction_velocity_ms, :u_star, :ustar,
        :u_star_5m, :ustar_5m, :u_star_m_s, :ustar_m_s, :u_star_tower, :ustar_tower
    ),
)

const EXPECTED_UNITS = Dict(
    :z => "m",
    :u => "m s^-1",
    :v => "m s^-1",
    :theta => "K",
    :q => "kg kg^-1",
    :u_star => "m s^-1",
)

const SURFACE_FLUX_ALIASES = (:hs, :h, :shf, :sensible_heat_flux, :h_sensible)
const THETA_REF_ALIASES = (:theta, :theta_k, :potential_temp, :potential_temperature, :tair, :temperature)
const UW_ALIASES = (:uw, :cov_uw, :u_w_, :uw_flux, :momentum_flux_u)
const VW_ALIASES = (:vw, :cov_vw, :v_w_, :vw_flux, :momentum_flux_v)

function _first_present_header(hindex::Dict{Symbol,Int}, aliases::Tuple)
    for a in aliases
        haskey(hindex, a) && return a
    end
    return nothing
end

function _derive_ustar_from_covariances(uw::Vector{Float64}, vw::Vector{Float64})
    n = min(length(uw), length(vw))
    out = fill(NaN, n)
    for i in 1:n
        if isfinite(uw[i]) && isfinite(vw[i])
            out[i] = (uw[i]^2 + vw[i]^2)^0.25
        end
    end
    return out
end

function _derive_obukhov_length(theta_ref::Vector{Float64}, ustar::Vector{Float64}, hs::Vector{Float64})
    n = min(length(theta_ref), length(ustar), length(hs))
    n == 0 && return Float64[]

    out = fill(NaN, n)
    kappa = 0.4
    rho_cp = 1200.0

    for i in 1:n
        θ = theta_ref[i]
        u = ustar[i]
        h = hs[i]
        if !(isfinite(θ) && isfinite(u) && isfinite(h)) || abs(h) <= 1e-12
            continue
        end
        out[i] = -θ * (u^3) * rho_cp / (kappa * h)
    end
    return out
end

function _to_f64(x)
    x isa Number && return Float64(x)
    parsed = tryparse(Float64, String(x))
    parsed === nothing && throw(ArgumentError("Encountered non-numeric value in tower dataset: $(x)"))
    return parsed
end

function _normalize_header(name)::Symbol
    return Symbol(replace(lowercase(String(name)), r"[^a-z0-9]+" => "_"))
end

function _normalize_column_selector(col)::Symbol
    col isa Symbol && return _normalize_header(col)
    return _normalize_header(String(col))
end

"""Resolve a raw CSV/DataFrame column name to its canonical variable symbol, if known."""
function resolve_header_alias(col)::Union{Symbol,Nothing}
    normalized = _normalize_header(col)
    for (canonical, aliases) in HEADER_ALIASES
        normalized in aliases && return canonical
    end
    return nothing
end

"""
    has_unit_suffix(col_name::Symbol) -> Bool

Returns true if the column header ends with a recognized unit suffix (e.g., _m, _ms, _k, _kgkg).
"""
function has_unit_suffix(col_name::Symbol)
    s = string(col_name)
    return any(endswith(s, suffix) for suffix in ["_m", "_ms", "_k", "_kgkg", "_wm2", "_pa", "_deg"])
end

"""Build an `ObservationTable` from a `DataFrame` using resolved canonical column mappings."""
function parse_observation_dataframe(df::DataFrame, mapped_cols::Dict{Symbol,Symbol}; strict_units::Bool=true, kwargs...)
    cols = Dict{Symbol,Vector{Float64}}()
    units = Dict{Symbol,String}()

    for (canonical, raw_col) in mapped_cols
        vals = [_to_f64(v) for v in df[!, raw_col]]
        strict_units && !all(isfinite, vals) && throw(ArgumentError("Column $(raw_col) contains non-finite values."))
        cols[canonical] = vals
        haskey(EXPECTED_UNITS, canonical) && (units[canonical] = EXPECTED_UNITS[canonical])
    end

    return ObservationTable(cols, units)
end

function _read_first_var(ds::NCDataset, names::Vector{Symbol})
    for n in names
        if haskey(ds, String(n))
            return vec(ds[String(n)][:])
        end
    end
    return nothing
end

function _read_required_var(ds::NCDataset, key::Union{Nothing,Symbol}, fallback_names::Vector{Symbol})
    if !(key === nothing) && haskey(ds, String(key))
        return vec(ds[String(key)][:])
    end
    res = _read_first_var(ds, fallback_names)
    res === nothing && throw(ArgumentError("None of the candidate NetCDF variables found: $(fallback_names)"))
    return res
end

"""
    read_tower_csv(filepath::AbstractString; strict_units::Bool = true, kwargs...)

Reads tower observation CSV data and validates required variable columns.
Throws an `ArgumentError` if required fields cannot be resolved from CSV headers.
"""
function read_tower_csv(filepath::AbstractString; strict_units::Bool = true, kwargs...)
    df = CSV.read(filepath, DataFrame)

    # Map raw CSV headers to canonical variable symbols
    mapped_cols = Dict{Symbol,Symbol}()
    for col in propertynames(df)
        if strict_units && !has_unit_suffix(col)
            continue # Skip bare unannotated headers under strict_units mode
        end
        canonical = resolve_header_alias(col)
        if canonical !== nothing
            mapped_cols[canonical] = col
        end
    end

    # Validate presence of mandatory canonical variables
    required_vars = (:z, :u, :v, :theta)
    missing_vars = filter(v -> !haskey(mapped_cols, v), required_vars)

    if !isempty(missing_vars)
        throw(ArgumentError("CSV file '$filepath' missing required canonical columns with unit annotations: $(missing_vars)"))
    end

    return parse_observation_dataframe(df, mapped_cols; strict_units=strict_units, kwargs...)
end

"""Read tower observation files with automatic flux and u_star fallback resolution."""
function read_observation_data(
    path::AbstractString;
    z_col::Union{Symbol,String}=:z,
    u_col::Union{Symbol,String}=:u,
    v_col::Union{Symbol,String}=:v,
    temp_col::Union{Symbol,String}=:theta,
    q_col::Union{Nothing,Symbol,String}=nothing,
    ustar_col::Union{Nothing,Symbol,String}=nothing,
    hs_col::Union{Nothing,Symbol,String}=nothing,
    theta_ref_col::Union{Nothing,Symbol,String}=nothing,
    include_derived_obukhov::Bool=true,
    auto_surface_flux_aliases::Bool=true,
    delim::Char=(','),
    default_q::Float64=0.0,
    default_u_star::Float64=0.3,
)
    ext = lowercase(splitext(path)[2])

    if ext == ".csv"
        raw, hdr = readdlm(path, delim, header=true)
        headers = [_normalize_header(h) for h in vec(hdr)]
        hindex = Dict(h => i for (i, h) in enumerate(headers))
        nrows = size(raw, 1)

        z_key = _normalize_column_selector(z_col)
        u_key = _normalize_column_selector(u_col)
        v_key = _normalize_column_selector(v_col)
        t_key = _normalize_column_selector(temp_col)
        q_key = q_col === nothing ? nothing : _normalize_column_selector(q_col)
        us_key = ustar_col === nothing ? nothing : _normalize_column_selector(ustar_col)
        hs_key = hs_col === nothing ? nothing : _normalize_column_selector(hs_col)
        θr_key = theta_ref_col === nothing ? nothing : _normalize_column_selector(theta_ref_col)

        if auto_surface_flux_aliases
            hs_key === nothing && (hs_key = _first_present_header(hindex, SURFACE_FLUX_ALIASES))
            θr_key === nothing && (θr_key = _first_present_header(hindex, THETA_REF_ALIASES))
            us_key === nothing && (us_key = _first_present_header(hindex, HEADER_ALIASES[:u_star]))
        end

        function col_from_raw(key::Symbol)
            j = hindex[key]
            return [_to_f64(raw[i, j]) for i in 1:nrows]
        end

        z = col_from_raw(z_key)
        u = col_from_raw(u_key)
        v = col_from_raw(v_key)
        theta = col_from_raw(t_key)
        q = q_key === nothing ? fill(default_q, nrows) : col_from_raw(q_key)

        # u_star resolution sequence: direct match -> kinematic covariance derivation -> fallback
        uw_key = _first_present_header(hindex, UW_ALIASES)
        vw_key = _first_present_header(hindex, VW_ALIASES)

        if us_key !== nothing
            ustar = col_from_raw(us_key)
        elseif uw_key !== nothing && vw_key !== nothing
            @info "Deriving u_star from covariance fluxes ($(uw_key), $(vw_key)) for $(path)"
            ustar = _derive_ustar_from_covariances(col_from_raw(uw_key), col_from_raw(vw_key))
        else
            @warn "u_star column missing in $(path). Using default fallback: $(default_u_star) m s^-1"
            ustar = fill(default_u_star, nrows)
        end

        hs = hs_key === nothing ? Float64[] : col_from_raw(hs_key)
        θr = θr_key === nothing ? theta : col_from_raw(θr_key)

        cols = Dict(:z => z, :u => u, :v => v, :theta => theta, :q => q, :u_star => ustar)
        if !isempty(hs)
            cols[:sensible_heat_flux] = hs
            include_derived_obukhov && (cols[:L_obukhov] = _derive_obukhov_length(θr, ustar, hs))
        end

        units = Dict(req => EXPECTED_UNITS[req] for req in REQUIRED_TOWER_COLUMNS)
        haskey(cols, :sensible_heat_flux) && (units[:sensible_heat_flux] = "W m^-2")
        haskey(cols, :L_obukhov) && (units[:L_obukhov] = "m")
        return ObservationTable(cols, units)

    elseif ext == ".nc"
        ds = NCDataset(path)
        try
            z = _read_required_var(ds, Symbol(z_col), [:z, :height, :zf, :zh])
            u = _read_required_var(ds, Symbol(u_col), [:u, :uwind, :u_component])
            v = _read_required_var(ds, Symbol(v_col), [:v, :vwind, :v_component])
            theta = _read_required_var(ds, Symbol(temp_col), [:theta, :potential_temperature, :th])
            q = q_col === nothing ? fill(default_q, length(z)) : _read_required_var(ds, Symbol(q_col), [:q, :specific_humidity])

            ustar_raw = _read_first_var(ds, [:u_star, :ustar, :friction_velocity, :u_star_5m, :ustar_tower, :u_star_ms])
            if ustar_raw !== nothing
                ustar = Float64.(ustar_raw)
            else
                uw = _read_first_var(ds, [:uw, :cov_uw, :u_w_])
                vw = _read_first_var(ds, [:vw, :cov_vw, :v_w_])
                if uw !== nothing && vw !== nothing
                    @info "Deriving u_star from NetCDF covariance fluxes (uw, vw) for $(path)"
                    ustar = _derive_ustar_from_covariances(Float64.(uw), Float64.(vw))
                else
                    @warn "u_star missing in NetCDF $(path). Using default fallback: $(default_u_star) m s^-1"
                    ustar = fill(default_u_star, length(z))
                end
            end

            n = minimum((length(z), length(u), length(v), length(theta), length(q), length(ustar)))
            cols = Dict(
                :z => Float64.(z[1:n]), :u => Float64.(u[1:n]), :v => Float64.(v[1:n]),
                :theta => Float64.(theta[1:n]), :q => Float64.(q[1:n]), :u_star => Float64.(ustar[1:n])
            )

            units = Dict(req => EXPECTED_UNITS[req] for req in REQUIRED_TOWER_COLUMNS)
            return ObservationTable(cols, units)
        finally
            close(ds)
        end
    else
        throw(ArgumentError("Unsupported extension $(ext). Use .csv or .nc."))
    end
end

"""Wrapper to read NetCDF tower observations into an ObservationTable."""
function read_tower_netcdf(path::AbstractString; kwargs...)
    return read_observation_data(path; kwargs...)
end

"""
    resolve_sibling_data_dir(target_dir::AbstractString = "data"; must_exist::Bool = false) -> String

Locates a sibling data directory relative to the package root or working directory.
If `must_exist` is true and no matching directory is found, throws an `ArgumentError`.
"""
function resolve_sibling_data_dir(target_dir::AbstractString = "data"; must_exist::Bool = false)
    pkg_root = normpath(joinpath(@__DIR__, "..", ".."))
    candidate = joinpath(pkg_root, target_dir)
    if isdir(candidate)
        return candidate
    end
    cwd_candidate = joinpath(pwd(), target_dir)
    if isdir(cwd_candidate)
        return cwd_candidate
    end
    if must_exist
        throw(ArgumentError("Directory '$target_dir' not found relative to package root or current working directory."))
    end
    return candidate
end

"""
    find_data_files(dir::AbstractString, ext::AbstractString = ".csv") -> Vector{String}

Recursively finds all data files matching the specified extension within a directory.
"""
function find_data_files(dir::AbstractString, ext::AbstractString = ".csv")
    !isdir(dir) && return String[]
    matches = String[]
    for (root, _, files) in walkdir(dir)
        for file in files
            if endswith(file, ext)
                push!(matches, joinpath(root, file))
            end
        end
    end
    return sort!(matches)
end