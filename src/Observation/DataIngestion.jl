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

function _normalize_header(name)::Symbol
    return Symbol(replace(lowercase(String(name)), r"[^a-z0-9]+" => "_"))
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
