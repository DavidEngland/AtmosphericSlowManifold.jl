using DelimitedFiles

struct ObservationTable
    columns::Dict{Symbol, Vector{Float64}}
end

function _normalize_header(name)
    return Symbol(replace(lowercase(String(name)), r"[^a-z0-9]+" => "_"))
end

"""Read tower profile CSV into a normalized in-memory observation table."""
function read_tower_csv(path::AbstractString; delim::Char = ',')
    raw, hdr = readdlm(path, delim, header = true)
    headers = [_normalize_header(h) for h in hdr]

    cols = Dict{Symbol, Vector{Float64}}()
    for (j, h) in enumerate(headers)
        cols[h] = [Float64(raw[i, j]) for i in 1:size(raw, 1)]
    end

    return ObservationTable(cols)
end

"""NetCDF support placeholder for phase-2 ingestion."""
function read_tower_netcdf(path::AbstractString)
    throw(ArgumentError("NetCDF ingestion is not implemented yet: $(path)"))
end
