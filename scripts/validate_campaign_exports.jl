using CSV
using DataFrames
using JSON3

const ROOT = normpath(joinpath(pwd(), "reports", "generated", "campaign_exports"))
const CSV_DIR = joinpath(ROOT, "csv")
const JSON_DIR = joinpath(ROOT, "json")
const NC_DIR = joinpath(ROOT, "netcdf")
const FIG_DIR = joinpath(ROOT, "figures")
const TABLE_DIR = joinpath(ROOT, "tables")

const CAMPAIGNS = ["cases_99", "floss", "bllast", "sheba"]

normalize_campaign(s::AbstractString) = lowercase(replace(s, r"[^A-Za-z0-9]+" => "_"))

function check(cond::Bool, msg::String)
    if cond
        println("[ok]   ", msg)
        return true
    else
        println("[fail] ", msg)
        return false
    end
end

function require_file(path::String)
    return check(isfile(path), "file exists: $(path)")
end

function require_dir(path::String)
    return check(isdir(path), "directory exists: $(path)")
end

function main()
    ok = true

    ok &= require_dir(ROOT)
    ok &= require_dir(CSV_DIR)
    ok &= require_dir(JSON_DIR)
    ok &= require_dir(NC_DIR)
    ok &= require_dir(FIG_DIR)
    ok &= require_dir(TABLE_DIR)

    summary_csv = joinpath(TABLE_DIR, "campaign_summary.csv")
    overview_csv = joinpath(TABLE_DIR, "campaign_overview.csv")
    summary_md = joinpath(TABLE_DIR, "campaign_summary.md")
    summary_tex = joinpath(TABLE_DIR, "campaign_summary.tex")

    ok &= require_file(summary_csv)
    ok &= require_file(overview_csv)
    ok &= require_file(summary_md)
    ok &= require_file(summary_tex)

    for c in CAMPAIGNS
        ok &= require_file(joinpath(CSV_DIR, "$(c)_raw.csv"))
        ok &= require_file(joinpath(CSV_DIR, "$(c)_stats.csv"))
        ok &= require_file(joinpath(JSON_DIR, "$(c)_model_and_diagnostics.json"))
        ok &= require_file(joinpath(FIG_DIR, "$(c)_metrics.png"))
    end

    # Validate JSON model/diagnostic schema for physical similarity parameters.
    for c in CAMPAIGNS
        json_path = joinpath(JSON_DIR, "$(c)_model_and_diagnostics.json")
        if !isfile(json_path)
            ok &= check(false, "JSON exists for schema validation: $(json_path)")
            continue
        end

        payload = JSON3.read(read(json_path, String))
        diag = payload["diagnostics"]

        has_similarity = haskey(diag, "similarity_parameters")
        ok &= check(has_similarity, "JSON diagnostics include similarity_parameters for $(c)")
        if has_similarity
            sim = diag["similarity_parameters"]
            has_zeta = haskey(sim, "zeta") && haskey(sim["zeta"], "n")
            has_phi = haskey(sim, "phi_obs") && haskey(sim["phi_obs"], "n")
            ok &= check(has_zeta, "similarity_parameters has zeta statistics for $(c)")
            ok &= check(has_phi, "similarity_parameters has phi_obs statistics for $(c)")
            if has_zeta
                ok &= check(Int(sim["zeta"]["n"]) > 0, "zeta finite record count N>0 for $(c)")
            end
            if has_phi
                ok &= check(Int(sim["phi_obs"]["n"]) > 0, "phi_obs finite record count N>0 for $(c)")
            end
        end

        has_obukhov = haskey(diag, "obukhov_scaling") && haskey(diag["obukhov_scaling"], "n")
        ok &= check(has_obukhov, "JSON diagnostics include obukhov_scaling.n for $(c)")
        if has_obukhov
            ok &= check(Int(diag["obukhov_scaling"]["n"]) > 0, "obukhov scaling finite record count N>0 for $(c)")
        end

        if haskey(payload, "terms")
            term_names = String[]
            for t in payload["terms"]
                haskey(t, "name") && push!(term_names, String(t["name"]))
            end
            has_abstract_a = any(startswith(n, "a_") for n in term_names)
            ok &= check(!has_abstract_a, "model terms exclude abstract a_* coefficients for $(c)")
        end
    end

    ok &= require_file(joinpath(NC_DIR, "cases_99_ri_profile.nc"))
    ok &= require_file(joinpath(NC_DIR, "floss_ri_profile.nc"))
    ok &= require_file(joinpath(NC_DIR, "bllast_ri_profile.nc"))
    ok &= require_file(joinpath(NC_DIR, "sheba_stability.nc"))

    # Heatmap is optional for SHEBA, required for trajectory campaigns.
    ok &= require_file(joinpath(FIG_DIR, "cases_99_ri_heatmap.png"))
    ok &= require_file(joinpath(FIG_DIR, "floss_ri_heatmap.png"))
    ok &= require_file(joinpath(FIG_DIR, "bllast_ri_heatmap.png"))
    ok &= require_file(joinpath(FIG_DIR, "campaign_overview.png"))

    # Validate overview schema.
    if isfile(overview_csv)
        overview = CSV.read(overview_csv, DataFrame)
        required_cols = [:campaign, :observations, :height_levels, :mean_wind_speed, :mean_ri]
        present_cols = Set(Symbol.(names(overview)))
        ok &= check(all(c -> c in present_cols, required_cols), "overview CSV has required columns")

        present = Set(normalize_campaign.(String.(overview.campaign)))
        expected = Set(CAMPAIGNS)
        ok &= check(present == expected, "overview CSV contains expected campaigns")

        ok &= check(all(overview.observations .> 0), "overview observations are positive")
        ok &= check(all(overview.height_levels .>= 0), "overview height levels are non-negative")
    end

    # Validate markdown report sections.
    if isfile(summary_md)
        txt = read(summary_md, String)
        ok &= check(occursin("## Campaign Analysis Summary", txt), "summary markdown contains analysis section")
        ok &= check(occursin("## Output Artifact Manifest", txt), "summary markdown contains artifact manifest section")
        ok &= check(occursin("Overview figure:", txt), "summary markdown references overview figure")
    end

    if ok
        println("Campaign artifact validation passed.")
        return
    end

    println("Campaign artifact validation failed.")
    exit(1)
end

main()
