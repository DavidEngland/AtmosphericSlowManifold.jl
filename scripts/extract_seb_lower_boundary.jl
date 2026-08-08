#!/usr/bin/env julia

using CSV
using DataFrames
using Dates

const DEFAULT_SOURCE = normpath(
    @__DIR__,
    "..",
    "..",
    "SpectralBL-Analytics",
    "data",
    "sheba",
    "raw",
    "ncar_eol_dee002994255",
    "prof_file_all6_ed_hd.txt",
)
const DEFAULT_OUTPUT = normpath(
    @__DIR__,
    "..",
    "reports",
    "generated",
    "campaign_exports",
    "csv",
    "seb_lower_boundary.csv",
)
const SHEBA_EPOCH = DateTime(1997, 1, 1)
const MILLISECONDS_PER_DAY = 86_400_000

function sheba_timestamp(julian_day::Real)
    return SHEBA_EPOCH + Millisecond(round(Int, (julian_day - 1) * MILLISECONDS_PER_DAY))
end

function extract_sheba(source::AbstractString)
    isfile(source) || error("SHEBA source file not found: $source")

    raw = CSV.read(
        source,
        DataFrame;
        delim='\t',
        header=1,
        skipto=3,
        missingstring=["999", "9999"],
        silencewarnings=true,
    )

    required = [:JD, :Tsfc]
    missing_columns = setdiff(required, propertynames(raw))
    isempty(missing_columns) || error("Missing required columns: $(join(missing_columns, ", "))")

    valid = dropmissing(select(raw, required))
    filter!(row -> isfinite(row.JD) && isfinite(row.Tsfc), valid)

    count = nrow(valid)
    return DataFrame(
        campaign=fill("sheba", count),
        timestamp=sheba_timestamp.(valid.JD),
        G=Vector{Union{Missing, Float64}}(missing, count),
        Ts=Float64.(valid.Tsfc),
        SWC=Vector{Union{Missing, Float64}}(missing, count),
    )
end

function main(args)
    length(args) <= 2 || error("usage: extract_seb_lower_boundary.jl [SOURCE] [OUTPUT]")
    source = length(args) >= 1 ? args[1] : DEFAULT_SOURCE
    output = length(args) >= 2 ? args[2] : DEFAULT_OUTPUT

    extracted = extract_sheba(source)
    mkpath(dirname(output))
    CSV.write(output, extracted; dateformat="yyyy-mm-ddTHH:MM:SS.sss")
    println("Wrote $(nrow(extracted)) SHEBA surface-temperature rows to $output")
end

main(ARGS)