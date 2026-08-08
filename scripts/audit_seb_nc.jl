#!/usr/bin/env julia
# scripts/audit_seb_nc.jl — scan a directory tree for NetCDF files and report
# variables that look like surface energy budget (SEB) fields: ground/soil
# heat flux (G), soil/snow/ice temperature, skin temperature, net/longwave/
# shortwave radiation, soil moisture, snow depth.
#
# Julia equivalent of the filename/header sections of scripts/audit_seb.sh,
# using NCDatasets.jl instead of ncdump so no external netcdf-bin install
# is required.
#
# Usage:
#   julia --project=. scripts/audit_seb_nc.jl [root_dir]
#
# Defaults to the sibling SpectralBL-Analytics/data directory if no root is
# given and it exists next to this repo, matching audit_seb.sh's default.
#
# Add the dependency once, from the repo root:
#   julia --project=. -e 'using Pkg; Pkg.add("NCDatasets")'

using NCDatasets

# --- token patterns, same intent as the FNAME_PAT / COL_PAT in audit_seb.sh
const NAME_PATTERNS = [
    r"soil"i, r"snow"i, r"ice"i, r"thermistor"i, r"ground.?heat"i,
    r"skin"i, r"surface"i, r"radiation"i, r"moist"i, r"^g$"i, r"^g_"i,
    r"^ts$"i, r"^ts_"i, r"^swc"i, r"^shf$"i, r"tskin"i, r"tsfc"i,
    r"t_sfc"i, r"tice"i, r"tsnw"i, r"tsnow"i, r"qsub"i, r"^rn$"i,
    r"netrad"i, r"^lwd$"i, r"^lwu$"i, r"^swd$"i, r"^swu$"i,
]

matches_seb(name::AbstractString) = any(p -> occursin(p, name), NAME_PATTERNS)

function default_root()
    here = @__DIR__
    sibling = joinpath(here, "..", "..", "SpectralBL-Analytics", "data")
    return isdir(sibling) ? normpath(sibling) : pwd()
end

function find_nc_files(root::AbstractString)
    files = String[]
    for (dirpath, _, filenames) in walkdir(root)
        for fn in filenames
            if lowercase(splitext(fn)[2]) == ".nc"
                push!(files, joinpath(dirpath, fn))
            end
        end
    end
    return sort(files)
end

function scan_file(path::AbstractString)
    hits = Tuple{String,String,String,String}[]  # (varname, long_name, standard_name, units)
    try
        NCDataset(path) do ds
            for (vname, var) in ds
                vstr = String(vname)
                long_name = haskey(var.attrib, "long_name") ? string(var.attrib["long_name"]) : ""
                standard_name = haskey(var.attrib, "standard_name") ? string(var.attrib["standard_name"]) : ""
                units = haskey(var.attrib, "units") ? string(var.attrib["units"]) : ""

                if any(matches_seb, (vstr, long_name, standard_name))
                    push!(hits, (vstr, long_name, standard_name, units))
                end
            end
        end
    catch e
        @warn "Could not open/read" path exception = e
    end
    return hits
end

function main()
    root = length(ARGS) >= 1 ? ARGS[1] : default_root()
    if !isdir(root)
        println(stderr, "Root directory not found: $root")
        exit(1)
    end
    println("Scanning for .nc files under: $root")
    nc_files = find_nc_files(root)
    println("Found $(length(nc_files)) NetCDF file(s).\n")

    println("############################################")
    println("# NetCDF variables matching SEB terms")
    println("############################################")

    any_hits = false
    for f in nc_files
        hits = scan_file(f)
        isempty(hits) && continue
        any_hits = true
        println("-- ", f, " --")
        for (vname, long_name, standard_name, units) in hits
            label = isempty(long_name) ? vname : "$vname ($long_name)"
            label = isempty(standard_name) ? label : "$label {standard_name=$standard_name}"
            label = isempty(units) ? label : "$label [$units]"
            println("    ", label)
        end
        println()
    end

    any_hits || println("  (no SEB-like variable names found)")

    println("Done. Cross-check hits against each campaign's README for exact")
    println("definitions, depths/heights, and QC flags before extraction.")
end

main()