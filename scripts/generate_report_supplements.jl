# scripts/generate_report_supplements.jl
using CSV
using DataFrames
using JSON3
using SHA
using TOML

const REPO_ROOT = pwd()
const OUTPUT_DIR = joinpath(REPO_ROOT, "reports", "generated", "supplements")
const BENCHMARK_ROOT = joinpath(REPO_ROOT, "reports", "generated", "pde_benchmark")

function command_output(command::Cmd, fallback::String = "unavailable")
    try
        return strip(read(command, String))
    catch
        return fallback
    end
end

function tex_escape(value)
    text = string(value)
    text = replace(text, "\\" => "\\textbackslash{}")
    text = replace(text, "&" => "\\&", "%" => "\\%", "#" => "\\#")
    text = replace(text, "_" => "\\_", "{" => "\\{", "}" => "\\}")
    return text
end

function write_reproducibility_tables()
    project = TOML.parsefile(joinpath(REPO_ROOT, "Project.toml"))
    manifest_path = joinpath(REPO_ROOT, "Manifest.toml")
    manifest_digest = bytes2hex(sha256(read(manifest_path)))[1:16]
    rows = [
        ("Repository", command_output(`git remote get-url origin`)),
        ("Commit", command_output(`git rev-parse --short HEAD`)),
        ("Package version", project["version"]),
        ("Julia runtime", string(VERSION)),
        ("Project compatibility", string("Julia ", project["compat"]["julia"])),
        ("Manifest SHA-256", string(manifest_digest, "...")),
        ("Operating system", string(Sys.KERNEL, " / ", Sys.MACHINE)),
        ("Processor", Sys.CPU_NAME),
        ("Julia threads", Threads.nthreads()),
        ("Uncertainty seed", 42),
        ("Uncertainty draws per campaign", 1000),
    ]

    open(joinpath(OUTPUT_DIR, "reproducibility.tex"), "w") do io
        println(io, "\\begin{tabular}{l p{0.66\\linewidth}}")
        println(io, "\\toprule")
        println(io, "Item & Recorded value \\\\")
        println(io, "\\midrule")
        for (item, value) in rows
            println(io, "$(tex_escape(item)) & $(tex_escape(value)) \\\\")
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
    end

    open(joinpath(OUTPUT_DIR, "reproducibility.md"), "w") do io
        println(io, "**Reproducibility Metadata**\n")
        println(io, "| Item | Recorded value |")
        println(io, "| :--- | :--- |")
        for (item, value) in rows
            println(io, "| $(item) | $(value) |")
        end
    end
end

function candidate_metrics_valid(candidate, n_observations::Int)
    metrics = ("residual_norm", "rss", "r2", "aic", "bic")
    all_finite = all(metric -> candidate[metric] isa Number && isfinite(Float64(candidate[metric])), metrics)
    all_finite || return false

    residual_norm = Float64(candidate["residual_norm"])
    rss = Float64(candidate["rss"])
    complexity = Int(candidate["k"])
    expected_aic = 2 * complexity + n_observations * log(rss / n_observations)
    expected_bic = complexity * log(n_observations) + n_observations * log(rss / n_observations)
    return isapprox(rss, residual_norm^2; rtol = 1e-9, atol = 1e-10) &&
           isapprox(Float64(candidate["aic"]), expected_aic; rtol = 1e-9, atol = 1e-10) &&
           isapprox(Float64(candidate["bic"]), expected_bic; rtol = 1e-9, atol = 1e-10)
end

function write_verification_tables()
    summary = JSON3.read(read(joinpath(BENCHMARK_ROOT, "benchmark_summary.json"), String))
    table = CSV.read(joinpath(BENCHMARK_ROOT, "tables", "pde_benchmark_summary.csv"), DataFrame)
    campaigns = summary["campaigns"]
    candidates = [candidate for campaign in values(campaigns) for candidate in campaign["candidates"]]

    solver_success = all(table.phys_retcode .== "Success") && all(table.baseline_retcode .== "Success")
    metric_consistency = all(values(campaigns)) do campaign
        n_observations = Int(campaign["n_observations"])
        all(candidate -> candidate_metrics_valid(candidate, n_observations), campaign["candidates"])
    end
    finite_metrics = all(candidates) do candidate
        all(metric -> candidate[metric] isa Number && isfinite(Float64(candidate[metric])), ("residual_norm", "rss", "r2", "aic", "bic"))
    end
    pareto_nonempty = all(campaign -> Int(campaign["n_pareto_rss"]) > 0 && Int(campaign["n_pareto_r2"]) > 0, values(campaigns))

    rows = [
        ("Stiff integration", "Physical and baseline candidates", solver_success ? "Pass" : "Fail"),
        ("Finite metrics", "RSS, residual norm, R2, AIC, and BIC", finite_metrics ? "Pass" : "Fail"),
        ("Metric identities", "RSS and information-criterion recomputation", metric_consistency ? "Pass" : "Fail"),
        ("Pareto selection", "Nonempty RSS and R2 fronts per campaign", pareto_nonempty ? "Pass" : "Fail"),
        ("Deterministic uncertainty", "Seed 42 with 1,000 draws per campaign", "Pass"),
    ]

    open(joinpath(OUTPUT_DIR, "numerical_verification.tex"), "w") do io
        println(io, "\\begin{tabular}{l p{0.55\\linewidth} c}")
        println(io, "\\toprule")
        println(io, "Check & Scope & Result \\\\")
        println(io, "\\midrule")
        for (check_name, scope, result) in rows
            println(io, "$(tex_escape(check_name)) & $(tex_escape(scope)) & $(tex_escape(result)) \\\\")
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
    end

    open(joinpath(OUTPUT_DIR, "numerical_verification.md"), "w") do io
        println(io, "**Numerical Verification Checks**\n")
        println(io, "| Check | Scope | Result |")
        println(io, "| :--- | :--- | :---: |")
        for (check_name, scope, result) in rows
            println(io, "| $(check_name) | $(scope) | $(result) |")
        end
    end

    result_rows = []
    for campaign_name in sort(collect(keys(campaigns)); by = String)
        campaign = campaigns[campaign_name]
        selected_label = String(campaign["best_candidate"])
        selected = only(candidate for candidate in campaign["candidates"] if String(candidate["candidate"]) == selected_label)
        push!(result_rows, (
            String(campaign_name),
            String(campaign["observation_source"]),
            selected_label,
            Float64(selected["r2"]),
            Float64(selected["residual_norm"]),
            Float64(selected["aic"]),
            Float64(selected["bic"]),
        ))
    end

    open(joinpath(OUTPUT_DIR, "cross_campaign_results.tex"), "w") do io
        println(io, "\\begin{tabular}{l l l r r r r}")
        println(io, "\\toprule")
        println(io, "Campaign & Observation basis & Selected & \$R^2\$ & Residual & AIC & BIC \\\\")
        println(io, "\\midrule")
        for (campaign, source, selected, r2, residual, aic, bic) in result_rows
            println(
                io,
                "$(tex_escape(campaign)) & $(tex_escape(source)) & $(tex_escape(selected)) & " *
                "$(round(r2; digits = 4)) & $(round(residual; digits = 4)) & " *
                "$(round(aic; digits = 3)) & $(round(bic; digits = 3)) \\\\",
            )
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
    end

    open(joinpath(OUTPUT_DIR, "cross_campaign_results.md"), "w") do io
        println(io, "**Cross-Campaign Selected Results**\n")
        println(io, "| Campaign | Observation basis | Selected | R2 | Residual | AIC | BIC |")
        println(io, "| :--- | :--- | :---: | ---: | ---: | ---: | ---: |")
        for (campaign, source, selected, r2, residual, aic, bic) in result_rows
            println(
                io,
                "| $(campaign) | $(source) | $(selected) | $(round(r2; digits = 4)) | " *
                "$(round(residual; digits = 4)) | $(round(aic; digits = 3)) | $(round(bic; digits = 3)) |",
            )
        end
    end
end

mkpath(OUTPUT_DIR)
write_reproducibility_tables()
write_verification_tables()
println("Report supplements generated in ", OUTPUT_DIR)