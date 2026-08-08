#!/usr/bin/env julia

using CSV
using DataFrames
using Dates
using NCDatasets

const REPO_ROOT = normpath(@__DIR__, "..")
const ANALYTICS_DATA = normpath(REPO_ROOT, "..", "SpectralBL-Analytics", "data")
const CASES99_DIR = joinpath(ANALYTICS_DATA, "cases99", "raw", "ncar_eol_dee0099881")
const SHEBA_FILE = joinpath(
    ANALYTICS_DATA,
    "sheba",
    "raw",
    "ncar_eol_dee002994255",
    "prof_file_all6_ed_hd.txt",
)
const OUT_PATH = joinpath(
    REPO_ROOT,
    "reports",
    "generated",
    "campaign_exports",
    "csv",
    "seb_lower_boundary.csv",
)
const SCHEMA = [
    :campaign,
    :timestamp,
    :station,
    :latitude,
    :longitude,
    :altitude_m,
    :level,
    :depth_cm,
    :G,
    :Ts,
    :SWC,
]
const SHEBA_EPOCH = DateTime(1997, 1, 1)
const MILLISECONDS_PER_DAY = 86_400_000

clean(::Missing) = missing
clean(value::Real) = !isfinite(value) || abs(value) >= 900 ? missing : Float64(value)

function clean_bounded(value, lower, upper)
    cleaned = clean(value)
    return ismissing(cleaned) || !(lower <= cleaned <= upper) ? missing : cleaned
end

function empty_output()
    return DataFrame(
        campaign=String[],
        timestamp=DateTime[],
        station=Union{Missing, String}[],
        latitude=Union{Missing, Float64}[],
        longitude=Union{Missing, Float64}[],
        altitude_m=Union{Missing, Float64}[],
        level=String[],
        depth_cm=Union{Missing, Float64}[],
        G=Union{Missing, Float64}[],
        Ts=Union{Missing, Float64}[],
        SWC=Union{Missing, Float64}[],
    )
end

function station_time_values(
    ds::NCDataset,
    name::AbstractString,
    station_count,
    time_count;
    bounds=(-Inf, Inf),
)
    haskey(ds, name) || return fill(missing, station_count, time_count)
    var = ds[name]
    dims = NCDatasets.dimnames(var)
    Set(dims) == Set(("station", "time")) || error("$name has unsupported dimensions $dims")

    values = clean_bounded.(Array(var[:, :]), bounds...)
    oriented = dims == ("station", "time") ? values : permutedims(values, (2, 1))
    size(oriented) == (station_count, time_count) ||
        error("$name has size $(size(oriented)); expected ($station_count, $time_count)")
    return oriented
end

function station_coordinate(ds::NCDataset, name::AbstractString, station_count)
    haskey(ds, name) || return fill(missing, station_count)
    values = clean.(Array(ds[name][:]))
    length(values) == station_count || error("$name does not match station dimension")
    return values
end

function extract_cases99(nc_dir::AbstractString)
    isdir(nc_dir) || error("CASES-99 directory not found: $nc_dir")
    files = sort(filter(path -> lowercase(splitext(path)[2]) == ".nc", readdir(nc_dir; join=true)))
    isempty(files) && error("No CASES-99 NetCDF files found in $nc_dir")
    rows = NamedTuple[]

    for path in files
        try
            NCDataset(path) do ds
                haskey(ds, "time") || error("missing time variable")
                times = Array(ds["time"][:])
                all(value -> ismissing(value) || value isa DateTime, times) ||
                    error("time did not decode to DateTime")

                haskey(ds, "latitude") || error("missing station latitude coordinate")
                station_count = length(ds["latitude"])
                time_count = length(times)
                latitude = station_coordinate(ds, "latitude", station_count)
                longitude = station_coordinate(ds, "longitude", station_count)
                altitude = station_coordinate(ds, "altitude", station_count)

                gsoil = station_time_values(ds, "Gsoil", station_count, time_count; bounds=(-500, 500))
                tsfc = station_time_values(ds, "Tsfc", station_count, time_count; bounds=(-80, 60))
                tsoil = station_time_values(ds, "Tsoil", station_count, time_count; bounds=(-80, 60))
                qsoil = station_time_values(ds, "Qsoil", station_count, time_count; bounds=(0, 1))

                for time_index in eachindex(times), station_index in 1:station_count
                    timestamp = times[time_index]
                    ismissing(timestamp) && continue
                    site = (
                        station="station_$(station_index)",
                        latitude=latitude[station_index],
                        longitude=longitude[station_index],
                        altitude_m=altitude[station_index],
                    )

                    surface_g = gsoil[station_index, time_index]
                    surface_t = tsfc[station_index, time_index]
                    if !ismissing(surface_g) || !ismissing(surface_t)
                        push!(rows, merge(
                            (campaign="cases_99", timestamp=timestamp),
                            site,
                            (level="surface", depth_cm=missing, G=surface_g, Ts=surface_t, SWC=missing),
                        ))
                    end

                    soil_t = tsoil[station_index, time_index]
                    soil_water = qsoil[station_index, time_index]
                    if !ismissing(soil_t) || !ismissing(soil_water)
                        push!(rows, merge(
                            (campaign="cases_99", timestamp=timestamp),
                            site,
                            (level="soil", depth_cm=missing, G=missing, Ts=soil_t, SWC=soil_water),
                        ))
                    end
                end
            end
        catch exception
            @warn "Skipping unreadable CASES-99 file" path exception
        end
    end

    return isempty(rows) ? empty_output() : DataFrame(rows)[:, SCHEMA]
end

function sheba_timestamp(julian_day::Real)
    return SHEBA_EPOCH + Millisecond(round(Int, (julian_day - 1) * MILLISECONDS_PER_DAY))
end

function extract_sheba(path::AbstractString)
    isfile(path) || error("SHEBA input file not found: $path")
    raw = CSV.read(
        path,
        DataFrame;
        delim='\t',
        header=1,
        skipto=3,
        missingstring=["999", "9999"],
        silencewarnings=true,
    )
    required = [:JD, :lat, :lon, :Tsfc]
    missing_columns = setdiff(required, propertynames(raw))
    isempty(missing_columns) || error("SHEBA is missing columns: $(join(missing_columns, ", "))")

    rows = NamedTuple[]
    for row in eachrow(raw)
        any(ismissing, (row.JD, row.Tsfc)) && continue
        surface_t = clean_bounded(row.Tsfc, -80, 40)
        ismissing(surface_t) && continue
        push!(rows, (
            campaign="sheba",
            timestamp=sheba_timestamp(row.JD),
            station="tower",
            latitude=clean(row.lat),
            longitude=clean(row.lon),
            altitude_m=missing,
            level="surface",
            depth_cm=missing,
            G=missing,
            Ts=surface_t,
            SWC=missing,
        ))
    end

    return isempty(rows) ? empty_output() : DataFrame(rows)[:, SCHEMA]
end

function main()
    println("Extracting CASES-99 station-resolved surface and soil fields...")
    cases99 = extract_cases99(CASES99_DIR)
    println("Extracting SHEBA surface temperature...")
    sheba = extract_sheba(SHEBA_FILE)

    combined = vcat(cases99, sheba)
    sort!(combined, [:campaign, :timestamp, :station, :level])
    mkpath(dirname(OUT_PATH))
    CSV.write(OUT_PATH, combined; dateformat="yyyy-mm-ddTHH:MM:SS.sss")

    println("Wrote $(nrow(combined)) rows to $OUT_PATH")
    println("  CASES-99 surface: ", count(row -> row.campaign == "cases_99" && row.level == "surface", eachrow(combined)))
    println("  CASES-99 soil:    ", count(row -> row.campaign == "cases_99" && row.level == "soil", eachrow(combined)))
    println("  SHEBA surface:    ", count(row -> row.campaign == "sheba", eachrow(combined)))
end

main()
