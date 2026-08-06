using DelimitedFiles
using NCDatasets

struct ObservationTable
    columns::Dict{Symbol, Vector{Float64}}
    units::Dict{Symbol, String}
end

const REQUIRED_TOWER_COLUMNS = [:z, :u, :v, :theta, :q, :u_star]

const HEADER_ALIASES = Dict(
    :z => (:z_m, :height_m, :z),
    :u => (:u_ms, :u_m_s, :u),
    :v => (:v_ms, :v_m_s, :v),
    :theta => (:theta_k, :potential_temperature_k, :theta),
    :q => (:q_kgkg, :specific_humidity_kgkg, :q),
    :u_star => (:u_star_ms, :ustar_ms, :friction_velocity_ms, :u_star, :ustar),
)

const EXPECTED_UNITS = Dict(
    :z => "m",
    :u => "m s^-1",
    :v => "m s^-1",
    :theta => "K",
    :q => "kg kg^-1",
    :u_star => "m s^-1",
)

const SURFACE_FLUX_ALIASES = (
    :hs,
    :h,
    :shf,
    :sensible_heat_flux,
    :h_sensible,
)

const THETA_REF_ALIASES = (
    :theta,
    :theta_k,
    :potential_temp,
    :potential_temperature,
    :tair,
    :temperature,
)

function _first_present_header(hindex::Dict{Symbol, Int}, aliases::Tuple)
    for a in aliases
        haskey(hindex, a) && return a
    end
    return nothing
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
        if !(isfinite(θ) && isfinite(u) && isfinite(h))
            continue
        end
        abs(h) <= 1e-12 && continue
        out[i] = -θ * (u^3) * rho_cp / (kappa * h)
    end
    return out
end

function resolve_sibling_data_dir(
    ;
    package_module::Module = AtmosphericSlowManifold,
    sibling_project::AbstractString = "SpectralBL-Analytics",
    data_subdir::AbstractString = "data",
    must_exist::Bool = true,
)
    data_dir = normpath(joinpath(pkgdir(package_module), "..", sibling_project, data_subdir))
    if must_exist && !isdir(data_dir)
        throw(ArgumentError("Directory does not exist: $(data_dir)"))
    end
    return data_dir
end

function find_data_files(
    root_dir::AbstractString;
    extensions::Vector{String} = [".csv", ".nc", ".h5"],
    recursive::Bool = true,
)
    isdir(root_dir) || throw(ArgumentError("Root directory does not exist: $(root_dir)"))
    extset = Set(lowercase.(extensions))
    found = String[]

    if recursive
        for (root, _, files) in walkdir(root_dir)
            for file in files
                ext = lowercase(splitext(file)[2])
                ext in extset || continue
                push!(found, joinpath(root, file))
            end
        end
    else
        for file in readdir(root_dir)
            fpath = joinpath(root_dir, file)
            isfile(fpath) || continue
            ext = lowercase(splitext(file)[2])
            ext in extset || continue
            push!(found, fpath)
        end
    end

    sort!(found)
    return found
end

function _normalize_header(name)::Symbol
    return Symbol(replace(lowercase(String(name)), r"[^a-z0-9]+" => "_"))
end

function _normalize_column_selector(col)::Symbol
    if col isa Symbol
        return _normalize_header(col)
    end
    return _normalize_header(String(col))
end

function _canonical_header_map(headers::Vector{Symbol})
    out = Dict{Symbol, Symbol}()
    for req in REQUIRED_TOWER_COLUMNS
        aliases = HEADER_ALIASES[req]
        idx = findfirst(h -> h in aliases, headers)
        idx === nothing && throw(ArgumentError("Missing required tower column for $(req). Allowed headers: $(collect(aliases))"))
        out[req] = headers[idx]
    end
    return out
end

function _strict_units_from_headers(cmap::Dict{Symbol, Symbol})
    units = Dict{Symbol, String}()
    for req in REQUIRED_TOWER_COLUMNS
        hdr = String(cmap[req])
        if endswith(hdr, "_m")
            units[req] = "m"
        elseif endswith(hdr, "_k")
            units[req] = "K"
        elseif endswith(hdr, "_kgkg")
            units[req] = "kg kg^-1"
        elseif endswith(hdr, "_ms") || endswith(hdr, "_m_s")
            units[req] = "m s^-1"
        else
            throw(ArgumentError("Column $(hdr) is missing strict unit suffix for $(req). Use headers like z_m, u_ms, theta_k, q_kgkg."))
        end
    end
    return units
end

function _assert_units(units::Dict{Symbol, String})
    for req in REQUIRED_TOWER_COLUMNS
        got = units[req]
        exp = EXPECTED_UNITS[req]
        got == exp || throw(ArgumentError("Unit mismatch for $(req): expected $(exp), got $(got)"))
    end
end

function _to_f64(x)
    if x isa Number
        return Float64(x)
    end
    parsed = tryparse(Float64, String(x))
    parsed === nothing && throw(ArgumentError("Encountered non-numeric value in tower dataset: $(x)"))
    return parsed
end

function _normalize_table(raw, headers::Vector{Symbol}, cmap::Dict{Symbol, Symbol})
    hindex = Dict(h => i for (i, h) in enumerate(headers))
    nrows = size(raw, 1)
    cols = Dict{Symbol, Vector{Float64}}()

    for req in REQUIRED_TOWER_COLUMNS
        src = cmap[req]
        j = hindex[src]
        vals = Vector{Float64}(undef, nrows)
        for i in 1:nrows
            vals[i] = _to_f64(raw[i, j])
        end
        all(isfinite, vals) || throw(ArgumentError("Column $(src) contains non-finite values."))
        cols[req] = vals
    end

    return cols
end

"""Read tower profile CSV into a normalized in-memory observation table."""
function read_tower_csv(path::AbstractString; delim::Char = ',', strict_units::Bool = true)
    raw, hdr = readdlm(path, delim, header = true)
    headers = [_normalize_header(h) for h in vec(hdr)]
    cmap = _canonical_header_map(headers)
    units = strict_units ? _strict_units_from_headers(cmap) : Dict(req => EXPECTED_UNITS[req] for req in REQUIRED_TOWER_COLUMNS)
    _assert_units(units)
    cols = _normalize_table(raw, headers, cmap)
    return ObservationTable(cols, units)
end

function _read_first_var(ds::NCDataset, names::Vector{Symbol})
    for n in names
        if haskey(ds, String(n))
            return vec(ds[String(n)][:])
        end
    end
    throw(ArgumentError("None of the candidate NetCDF variables found: $(names)"))
end

function _read_required_var(ds::NCDataset, key::Union{Nothing, Symbol}, fallback_names::Vector{Symbol})
    if !(key === nothing) && haskey(ds, String(key))
        return vec(ds[String(key)][:])
    end
    return _read_first_var(ds, fallback_names)
end

"""Read tower/radiosonde/LES profile variables from NetCDF into normalized observation columns."""
function read_tower_netcdf(path::AbstractString)
    ds = NCDataset(path)
    try
        z = _read_first_var(ds, [:z, :height, :zf, :zh])
        u = _read_first_var(ds, [:u, :uwind, :u_component])
        v = _read_first_var(ds, [:v, :vwind, :v_component])
        theta = _read_first_var(ds, [:theta, :potential_temperature, :th])
        q = _read_first_var(ds, [:q, :specific_humidity])
        u_star = _read_first_var(ds, [:u_star, :ustar, :friction_velocity])

        n = minimum((length(z), length(u), length(v), length(theta), length(q), length(u_star)))
        n == 0 && throw(ArgumentError("No profile samples found in NetCDF file: $(path)"))

        cols = Dict(
            :z => Float64.(z[1:n]),
            :u => Float64.(u[1:n]),
            :v => Float64.(v[1:n]),
            :theta => Float64.(theta[1:n]),
            :q => Float64.(q[1:n]),
            :u_star => Float64.(u_star[1:n]),
        )

        for req in REQUIRED_TOWER_COLUMNS
            all(isfinite, cols[req]) || throw(ArgumentError("NetCDF variable $(req) contains non-finite values."))
        end

        units = Dict(req => EXPECTED_UNITS[req] for req in REQUIRED_TOWER_COLUMNS)
        return ObservationTable(cols, units)
    finally
        close(ds)
    end
end

"""
Read observation data with explicit column mapping.

Supports `.csv` and `.nc` files and normalizes output to `ObservationTable`
with canonical columns: `:z`, `:u`, `:v`, `:theta`, `:q`, `:u_star`.
"""
function read_observation_data(
    path::AbstractString;
    z_col::Union{Symbol, String} = :z,
    u_col::Union{Symbol, String} = :u,
    v_col::Union{Symbol, String} = :v,
    temp_col::Union{Symbol, String} = :theta,
    q_col::Union{Nothing, Symbol, String} = nothing,
    ustar_col::Union{Nothing, Symbol, String} = nothing,
    hs_col::Union{Nothing, Symbol, String} = nothing,
    theta_ref_col::Union{Nothing, Symbol, String} = nothing,
    include_derived_obukhov::Bool = true,
    auto_surface_flux_aliases::Bool = true,
    compute_obukhov::Union{Nothing, Bool} = nothing,
    surface_flux_aliases::Union{Nothing, Bool} = nothing,
    delim::Char = ',',
    default_q::Float64 = 0.0,
    default_u_star::Float64 = 0.3,
)
    if !(compute_obukhov === nothing)
        include_derived_obukhov = compute_obukhov
    end
    if !(surface_flux_aliases === nothing)
        auto_surface_flux_aliases = surface_flux_aliases
    end

    ext = lowercase(splitext(path)[2])

    if ext == ".csv"
        raw, hdr = readdlm(path, delim, header = true)
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

        for key in (z_key, u_key, v_key, t_key)
            haskey(hindex, key) || throw(ArgumentError("Missing required column $(key) in $(path)."))
        end

        function col_from_raw(key::Symbol)
            j = hindex[key]
            out = Vector{Float64}(undef, nrows)
            for i in 1:nrows
                out[i] = _to_f64(raw[i, j])
            end
            return out
        end

        z = col_from_raw(z_key)
        u = col_from_raw(u_key)
        v = col_from_raw(v_key)
        theta = col_from_raw(t_key)
        q = q_key === nothing ? fill(default_q, nrows) : col_from_raw(q_key)
        ustar = us_key === nothing ? fill(default_u_star, nrows) : col_from_raw(us_key)
        hs = hs_key === nothing ? Float64[] : col_from_raw(hs_key)
        θr = θr_key === nothing ? theta : col_from_raw(θr_key)

        cols = Dict(
            :z => z,
            :u => u,
            :v => v,
            :theta => theta,
            :q => q,
            :u_star => ustar,
        )
        if !isempty(hs)
            cols[:sensible_heat_flux] = hs
        end
        if include_derived_obukhov && !isempty(hs)
            cols[:L_obukhov] = _derive_obukhov_length(θr, ustar, hs)
        end
        for req in REQUIRED_TOWER_COLUMNS
            all(isfinite, cols[req]) || throw(ArgumentError("Column $(req) contains non-finite values in $(path)."))
        end
        units = Dict(req => EXPECTED_UNITS[req] for req in REQUIRED_TOWER_COLUMNS)
        if haskey(cols, :sensible_heat_flux)
            units[:sensible_heat_flux] = "W m^-2"
        end
        if haskey(cols, :L_obukhov)
            units[:L_obukhov] = "m"
        end
        return ObservationTable(cols, units)
    elseif ext == ".nc"
        ds = NCDataset(path)
        try
            z = _read_required_var(ds, Symbol(z_col), [:z, :height, :zf, :zh])
            u = _read_required_var(ds, Symbol(u_col), [:u, :uwind, :u_component])
            v = _read_required_var(ds, Symbol(v_col), [:v, :vwind, :v_component])
            theta = _read_required_var(ds, Symbol(temp_col), [:theta, :potential_temperature, :th])
            q = q_col === nothing ? fill(default_q, length(z)) : _read_required_var(ds, Symbol(q_col), [:q, :specific_humidity])
            ustar = ustar_col === nothing ? fill(default_u_star, length(z)) : _read_required_var(ds, Symbol(ustar_col), [:u_star, :ustar, :friction_velocity])
            hs = hs_col === nothing ? Float64[] : _read_required_var(ds, Symbol(hs_col), [:hs, :h, :shf, :sensible_heat_flux])
            θr = theta_ref_col === nothing ? theta : _read_required_var(ds, Symbol(theta_ref_col), [:theta, :potential_temperature, :th, :temperature])

            if auto_surface_flux_aliases && isempty(hs)
                try
                    hs = _read_first_var(ds, [:hs, :h, :shf, :sensible_heat_flux])
                catch
                end
            end

            n = minimum((length(z), length(u), length(v), length(theta), length(q), length(ustar)))
            n == 0 && throw(ArgumentError("No profile samples found in NetCDF file: $(path)"))

            cols = Dict(
                :z => Float64.(z[1:n]),
                :u => Float64.(u[1:n]),
                :v => Float64.(v[1:n]),
                :theta => Float64.(theta[1:n]),
                :q => Float64.(q[1:n]),
                :u_star => Float64.(ustar[1:n]),
            )
            if !isempty(hs)
                m = min(n, length(hs))
                cols[:sensible_heat_flux] = Float64.(hs[1:m])
                if include_derived_obukhov
                    mθ = min(n, length(θr), length(ustar), length(hs))
                    cols[:L_obukhov] = _derive_obukhov_length(
                        Float64.(θr[1:mθ]),
                        Float64.(ustar[1:mθ]),
                        Float64.(hs[1:mθ]),
                    )
                end
            end
            for req in REQUIRED_TOWER_COLUMNS
                all(isfinite, cols[req]) || throw(ArgumentError("NetCDF variable $(req) contains non-finite values."))
            end
            units = Dict(req => EXPECTED_UNITS[req] for req in REQUIRED_TOWER_COLUMNS)
            if haskey(cols, :sensible_heat_flux)
                units[:sensible_heat_flux] = "W m^-2"
            end
            if haskey(cols, :L_obukhov)
                units[:L_obukhov] = "m"
            end
            return ObservationTable(cols, units)
        finally
            close(ds)
        end
    else
        throw(ArgumentError("Unsupported observation file extension $(ext). Use .csv or .nc."))
    end
end
