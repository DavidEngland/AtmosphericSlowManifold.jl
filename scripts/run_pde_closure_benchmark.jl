# scripts/run_pde_closure_benchmark.jl
using AtmosphericSlowManifold
using CSV
using DataFrames
using Dates
using JSON3
using Plots
using Plots.Measures
using Statistics
using DifferentialEquations

const INPUT_ROOT = joinpath(pwd(), "reports", "generated", "campaign_exports")
const JSON_DIR = joinpath(INPUT_ROOT, "json")
const CSV_DIR = joinpath(INPUT_ROOT, "csv")

const OUTPUT_ROOT = joinpath(pwd(), "reports", "generated", "pde_benchmark")
const SUMMARY_PATH = joinpath(OUTPUT_ROOT, "benchmark_summary.json")
const FIG_DIR = joinpath(OUTPUT_ROOT, "figures")
const TABLE_DIR = joinpath(OUTPUT_ROOT, "tables")

const N_MODES_BENCHMARK = 12
const RELTOL_BENCHMARK = 1e-6
const ABSTOL_BENCHMARK = 1e-8
const SAVEAT_BENCHMARK = 300.0
const DT_INITIAL_BENCHMARK = 1.0
const MAXITERS_BENCHMARK = 10^7
const NONLINEAR_SCALE_BENCHMARK = 0.1

mkpath(FIG_DIR)
mkpath(TABLE_DIR)

const CAMPAIGNS = Dict(
    "SHEBA" => (
        json = joinpath(JSON_DIR, "sheba_model_and_diagnostics.json"),
        raw = joinpath(CSV_DIR, "sheba_raw.csv"),
    ),
    "CASES-99" => (
        json = joinpath(JSON_DIR, "cases_99_model_and_diagnostics.json"),
        raw = joinpath(CSV_DIR, "cases_99_raw.csv"),
    ),
)

function to_f64(x)
    if x === missing || x === nothing
        return NaN
    elseif x isa Number
        return Float64(x)
    end
    v = tryparse(Float64, strip(String(x)))
    return v === nothing ? NaN : v
end

function finite_mean(v)
    vf = filter(isfinite, v)
    return isempty(vf) ? NaN : mean(vf)
end

function campaign_mean_wind(raw_csv::String)
    df = CSV.read(raw_csv, DataFrame)
    rename!(df, Symbol.(names(df)))

    if hasproperty(df, :ws_lo) && hasproperty(df, :ws_hi)
        ws_lo = to_f64.(df[!, :ws_lo])
        ws_hi = to_f64.(df[!, :ws_hi])
        n = min(length(ws_lo), length(ws_hi))
        return finite_mean((ws_lo[1:n] .+ ws_hi[1:n]) ./ 2)
    elseif hasproperty(df, :eta_1) && hasproperty(df, :eta_2)
        eta1 = to_f64.(df[!, :eta_1])
        eta2 = to_f64.(df[!, :eta_2])
        n = min(length(eta1), length(eta2))
        return finite_mean(hypot.(eta1[1:n], eta2[1:n]))
    elseif hasproperty(df, :ustar)
        return finite_mean(to_f64.(df[!, :ustar]))
    end

    return 1.0
end

function campaign_observed_modes(raw_csv::String, n_modes::Int, lambda::Float64)
    df = CSV.read(raw_csv, DataFrame)
    rename!(df, Symbol.(names(df)))
    modes = zeros(Float64, n_modes)

    modal_columns = Tuple{Int, Symbol}[]
    for name in propertynames(df)
        matched = match(r"^a_(\d+)$", String(name))
        matched === nothing || push!(modal_columns, (parse(Int, matched.captures[1]), name))
    end
    sort!(modal_columns; by = first)

    if !isempty(modal_columns)
        for (mode_index, name) in modal_columns
            mode_index >= n_modes && continue
            modes[mode_index + 1] = finite_mean(to_f64.(df[!, name]))
        end
        all(isfinite, modes) && return modes, "campaign_modal_coefficients"
    end

    if all(hasproperty(df, name) for name in (:ws_lo, :ws_hi))
        wind_lo = finite_mean(to_f64.(df[!, :ws_lo]))
        wind_hi = finite_mean(to_f64.(df[!, :ws_hi]))
        if isfinite(wind_lo) && isfinite(wind_hi)
            modes[1] = (wind_lo + wind_hi) / 2
            n_modes >= 2 && (modes[2] = (wind_hi - wind_lo) / (4 * lambda))
            return modes, "two_level_linear_profile"
        end
    end

    modes[1] = campaign_mean_wind(raw_csv)
    return modes, "campaign_mean_wind"
end

function rmse(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    length(a) == length(b) || throw(DimensionMismatch("Vectors must match for RMSE."))
    return sqrt(mean((a .- b) .^ 2))
end

function residual_norm(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    length(a) == length(b) || throw(DimensionMismatch("Vectors must match for residual norm."))
    return sqrt(sum(abs2, a .- b))
end

finite_or_missing(x::Real) = isfinite(x) ? Float64(x) : nothing

function finite_or_na_str(x::Real; digits::Int = 4)
    return isfinite(x) ? string(round(Float64(x); digits = digits)) : "NA"
end

function r2_score(pred::AbstractVector{<:Real}, truth::AbstractVector{<:Real})
    length(pred) == length(truth) || throw(DimensionMismatch("Vectors must match for R2."))
    tbar = mean(truth)
    rss = sum((pred .- truth) .^ 2)
    tss = sum((truth .- tbar) .^ 2)
    return tss <= 1e-12 ? (rss <= 1e-12 ? 1.0 : 0.0) : 1.0 - rss / tss
end

function write_text(path::String, content::String)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, content)
    end
    return path
end

function escape_tex(str::AbstractString)
    escaped = String(str)
    escaped = replace(escaped, "\\" => "\\textbackslash{}")
    escaped = replace(escaped, "_" => "\\_")
    escaped = replace(escaped, "%" => "\\%")
    escaped = replace(escaped, "&" => "\\&")
    escaped = replace(escaped, "#" => "\\#")
    escaped = replace(escaped, "{" => "\\{")
    escaped = replace(escaped, "}" => "\\}")
    escaped = replace(escaped, string('$') => string('\\', '$'))
    return escaped
end

function export_cross_campaign_summary_tex(
    filepath::String,
    summary_df::DataFrame,
    best_models::Dict{Symbol, DiscoveredModel{Float64}},
    r2_by_campaign::Dict{Symbol, Float64},
    residual_by_campaign::Dict{Symbol, Float64},
)
    open(filepath, "w") do io
        println(io, "\\begin{tabular}{l r r r}")
        println(io, "\\toprule")
        println(io, "Site & Active terms & " * "\$R^{2}\$" * " & Residual norm " * "\\\\")
        println(io, "\\midrule")

        for row in eachrow(summary_df)
            camp_sym = campaign_symbol(String(row.campaign))
            site_clean = escape_tex(string(camp_sym))
            k = get(best_models, camp_sym, nothing)
            active_terms = k === nothing ? 0 : k.sparsity_level
            r2_val = get(r2_by_campaign, camp_sym, NaN)
            resid_val = get(residual_by_campaign, camp_sym, NaN)

            if isfinite(r2_val)
                r2_str = string(round(r2_val; digits = 4))
            else
                r2_str = "NA"
            end
            if isfinite(resid_val)
                resid_str = string(round(resid_val; digits = 5))
            else
                resid_str = "NA"
            end

            println(io, "\\texttt{$(site_clean)} & $(active_terms) & $(r2_str) & $(resid_str) \\\\")
        end

        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
    end
    return filepath
end

function campaign_symbol(campaign::String)
    return Symbol(replace(lowercase(campaign), "-" => "_"))
end

function _site_slug(site_name::String)
    slug = lowercase(strip(site_name))
    slug = replace(slug, "-" => "_")
    slug = replace(slug, " " => "_")
    return slug
end

function model_from_profile(
    target::Symbol,
    mean_diffusivity::Float64,
    residual_norm::Float64,
)
    terms = OperatorTerm{Float64}[
        OperatorTerm{Float64}(mean_diffusivity, BasisOperator[BasisOperator(SpatialDerivative(:u, 2), 1.0)]),
    ]
    return DiscoveredModel{Float64}(target, terms, max(residual_norm, 1e-12), length(terms))
end

function make_pareto_plot(
    complexity::Vector{Int},
    score::Vector{Float64},
    labels::Vector{String},
    pareto_idx::Vector{Int},
    title_text::String,
    y_label::String,
    out_path::String,
)
    finite_idx = findall(isfinite, score)
    p = plot(
        xlabel = "Model complexity (k)",
        ylabel = y_label,
        title = title_text,
    )

    if isempty(finite_idx)
        annotate!(p, 0.5, 0.5, text("No finite candidate metrics", 10, :center))
        savefig(p, out_path)
        return out_path
    end

    scatter!(
        p,
        complexity[finite_idx],
        score[finite_idx],
        markersize = 7,
        label = "Candidates",
    )
    if !isempty(pareto_idx)
        pareto_finite = [i for i in pareto_idx if i in finite_idx]
        if !isempty(pareto_finite)
        scatter!(
            p,
            complexity[pareto_finite],
            score[pareto_finite],
            markersize = 9,
            markerstrokewidth = 2,
            markercolor = :red,
            label = "Pareto front",
        )
        end
    end
    for i in finite_idx
        annotate!(p, complexity[i], score[i], text(labels[i], 8, :left))
    end
    savefig(p, out_path)
    return out_path
end

function solve_with_fallback(pde_sys, closure, disc, tspan, u0)
    solver_chain = (
        (name = "Rodas5P", alg = Rodas5P(autodiff = true)),
        (name = "FBDF", alg = FBDF()),
        (name = "RadauIIA5", alg = RadauIIA5()),
    )

    failures = String[]
    for candidate in solver_chain
        try
            sol = solve_scm(
                pde_sys,
                closure,
                disc,
                tspan;
                solver = candidate.alg,
                reltol = RELTOL_BENCHMARK,
                abstol = ABSTOL_BENCHMARK,
                dt = DT_INITIAL_BENCHMARK,
                dtmax = 120.0,
                maxiters = MAXITERS_BENCHMARK,
                saveat = SAVEAT_BENCHMARK,
                u0 = u0,
                verbose = false,
            )

            if SciMLBase.successful_retcode(sol)
                return sol, candidate.name, true
            end
            push!(failures, "$(candidate.name): $(sol.retcode)")
        catch err
            push!(failures, "$(candidate.name): $(sprint(showerror, err))")
        end
    end

    @warn "All stiff solver candidates failed" failures
    return nothing, "failed", false
end

function run_case(campaign::String, payload)
    isfile(payload.json) || throw(ArgumentError("Missing JSON payload: $(payload.json)"))
    isfile(payload.raw) || throw(ArgumentError("Missing raw campaign CSV: $(payload.raw)"))

    closure = PhysicalSimilarityClosure(payload.json)
    baseline = PhysicalSimilarityClosure(
        phi_coeffs = [1.0],
        zeta_coeffs = [0.0, 1.0],
        karman = closure.karman,
        ustar = closure.ustar,
        L_obukhov = closure.L_obukhov,
        z_ref = closure.z_ref,
    )

    disc = SpectralBLGalerkin(
        n_modes = N_MODES_BENCHMARK,
        lambda = 0.75,
        H = 100.0,
        enable_nonlinear = true,
        nonlinear_scale = NONLINEAR_SCALE_BENCHMARK,
    )
    tspan = (0.0, 12.0 * 3600.0)

    observed, observation_source = campaign_observed_modes(payload.raw, disc.n_modes, disc.lambda)
    mean_wind = campaign_mean_wind(payload.raw)
    mean_wind = isfinite(mean_wind) ? mean_wind : observed[1]
    u0 = copy(observed)

    pde_phys = build_pde_system(closure; z_top = 100.0, t_end = tspan[2], coriolis = 0.0, v_geostrophic = 0.0, radiation = 0.0)
    pde_base = build_pde_system(baseline; z_top = 100.0, t_end = tspan[2], coriolis = 0.0, v_geostrophic = 0.0, radiation = 0.0)

    sol_phys, solver_phys, ok_phys = solve_with_fallback(pde_phys, closure, disc, tspan, u0)
    sol_base, solver_base, ok_base = solve_with_fallback(pde_base, baseline, disc, tspan, u0)

    final_phys = ok_phys ? Float64.(sol_phys.u[end]) : fill(NaN, disc.n_modes)
    final_base = ok_base ? Float64.(sol_base.u[end]) : fill(NaN, disc.n_modes)

    z_grid = collect(range(2.0, 100.0; length = disc.n_modes))
    K_phys = zeros(Float64, disc.n_modes)
    K_base = zeros(Float64, disc.n_modes)
    evaluate_diffusivity_profile!(K_phys, closure, z_grid)
    evaluate_diffusivity_profile!(K_base, baseline, z_grid)

    r2_phys = ok_phys ? r2_score(final_phys, observed) : -Inf
    r2_base = ok_base ? r2_score(final_base, observed) : -Inf

    baseline_model = model_from_profile(:u, mean(K_base), ok_base ? residual_norm(final_base, observed) : Inf)
    physical_model = model_from_profile(:u, mean(K_phys), ok_phys ? residual_norm(final_phys, observed) : Inf)

    candidates = [baseline_model, physical_model]
    closure_labels = ["baseline", "physical"]
    n_obs = length(u0)

    aic_vals = [model_aic(m, n_obs) for m in candidates]
    bic_vals = [model_bic(m, n_obs) for m in candidates]
    rss_vals = [m.residual_norm^2 for m in candidates]
    r2_vals = [r2_base, r2_phys]

    pareto_rss = compute_pareto_front(candidates; objective = :rss)
    pareto_r2 = compute_pareto_front(candidates; objective = :r2, r2_values = r2_vals)

    best_idx = if !isempty(pareto_rss.indices)
        pareto_rss.indices[1]
    else
        finite_idx = findall(isfinite, rss_vals)
        isempty(finite_idx) ? 1 : finite_idx[argmin(rss_vals[finite_idx])]
    end
    best_model = candidates[best_idx]
    best_label = closure_labels[best_idx]

    return (
        campaign = campaign,
        observation_source = observation_source,
        mean_wind = mean_wind,
        rmse_physical = ok_phys ? rmse(final_phys, observed) : Inf,
        rmse_baseline = ok_base ? rmse(final_base, observed) : Inf,
        r2_physical = r2_phys,
        r2_baseline = r2_base,
        max_km = maximum(K_phys),
        solver_physical = solver_phys,
        solver_baseline = solver_base,
        phys_retcode = ok_phys ? string(sol_phys.retcode) : "solve_failed",
        baseline_retcode = ok_base ? string(sol_base.retcode) : "solve_failed",
        final_phys = final_phys,
        final_base = final_base,
        observed = observed,
        candidate_models = candidates,
        closure_labels = closure_labels,
        aic_vals = aic_vals,
        bic_vals = bic_vals,
        rss_vals = rss_vals,
        r2_vals = r2_vals,
        pareto_rss = pareto_rss,
        pareto_r2 = pareto_r2,
        best_model = best_model,
        best_label = best_label,
        best_idx = best_idx,
    )
end

rows = NamedTuple[]
profiles = Dict{String, NamedTuple}()
best_models_by_campaign = Dict{Symbol, DiscoveredModel{Float64}}()
r2_best_by_campaign = Dict{Symbol, Float64}()
residual_best_by_campaign = Dict{Symbol, Float64}()
summary_campaigns = Dict{String, Any}()

for (campaign, payload) in CAMPAIGNS
    result = run_case(campaign, payload)
    camp_sym = campaign_symbol(campaign)

    push!(rows, (
        campaign = result.campaign,
        observation_source = result.observation_source,
        mean_wind = result.mean_wind,
        rmse_physical = result.rmse_physical,
        rmse_baseline = result.rmse_baseline,
        r2_physical = result.r2_physical,
        r2_baseline = result.r2_baseline,
        max_km = result.max_km,
        phys_retcode = result.phys_retcode,
        baseline_retcode = result.baseline_retcode,
        best_candidate = result.best_label,
        best_aic = result.aic_vals[result.best_idx],
        best_bic = result.bic_vals[result.best_idx],
    ))

    candidate_rows = Any[]
    for i in eachindex(result.candidate_models)
        push!(candidate_rows, Dict(
            "candidate" => result.closure_labels[i],
            "k" => result.candidate_models[i].sparsity_level,
            "residual_norm" => finite_or_missing(result.candidate_models[i].residual_norm),
            "rss" => finite_or_missing(result.rss_vals[i]),
            "r2" => finite_or_missing(result.r2_vals[i]),
            "aic" => finite_or_missing(result.aic_vals[i]),
            "bic" => finite_or_missing(result.bic_vals[i]),
            "pareto_rss" => in(i, result.pareto_rss.indices),
            "pareto_r2" => in(i, result.pareto_r2.indices),
        ))
    end

    eq_tex = to_latex(result.best_model)
    write_text(joinpath(TABLE_DIR, "$(camp_sym)_best_equation.tex"), string("\\[\n", eq_tex, "\n\\]\n"))
    write_text(
        joinpath(TABLE_DIR, "$(camp_sym)_best_terms.tex"),
        latex_term_table(
            result.best_model;
            r2 = result.r2_vals[result.best_idx],
            residual_norm = result.best_model.residual_norm,
        ),
    )

    complexities = [m.sparsity_level for m in result.candidate_models]
    make_pareto_plot(
        complexities,
        result.rss_vals,
        result.closure_labels,
        result.pareto_rss.indices,
        "$(campaign): Pareto (RSS vs k)",
        "RSS",
        joinpath(FIG_DIR, "$(camp_sym)_pareto_rss.png"),
    )
    make_pareto_plot(
        complexities,
        result.r2_vals,
        result.closure_labels,
        result.pareto_r2.indices,
        "$(campaign): Pareto (R2 vs k)",
        "R2",
        joinpath(FIG_DIR, "$(camp_sym)_pareto_r2.png"),
    )

    best_models_by_campaign[camp_sym] = result.best_model
    r2_best_by_campaign[camp_sym] = result.r2_vals[result.best_idx]
    residual_best_by_campaign[camp_sym] = result.best_model.residual_norm

    summary_campaigns[campaign] = Dict(
        "n_candidates" => length(result.candidate_models),
        "n_observations" => length(result.observed),
        "n_pareto_rss" => length(result.pareto_rss.indices),
        "n_pareto_r2" => length(result.pareto_r2.indices),
        "best_candidate" => result.best_label,
        "observation_source" => result.observation_source,
        "best_model" => Dict(
            "terms" => length(result.best_model.terms),
            "residual_norm" => finite_or_missing(result.best_model.residual_norm),
            "aic" => finite_or_missing(result.aic_vals[result.best_idx]),
            "bic" => finite_or_missing(result.bic_vals[result.best_idx]),
            "latex" => eq_tex,
        ),
        "candidates" => candidate_rows,
    )

    profiles[campaign] = result
end

summary_df = DataFrame(rows)
summary_path = joinpath(TABLE_DIR, "pde_benchmark_summary.csv")
CSV.write(summary_path, summary_df)

site_table_path = joinpath(TABLE_DIR, "cross_campaign_best_models.tex")
export_cross_campaign_summary_tex(
    site_table_path,
    summary_df,
    best_models_by_campaign,
    r2_best_by_campaign,
    residual_best_by_campaign,
)

artifact_root = joinpath(pwd(), "reports", "generated", "campaign_exports")
loso_discovery = function (_train_campaigns, X, y; kwargs...)
    coeffs = X \ y
    rss = sum(abs2, X * coeffs .- y)
    return (coefficients = coeffs, discovered_model = "linear_ols", train_rss = rss)
end

loso_summary = run_artifact_loso(
    artifact_root,
    loso_discovery;
    sites = ["SHEBA", "CASES-99", "FLOSS", "BLLAST"],
    feature_cols = [1],
    target_col = 2,
)

loso_table_path = joinpath(TABLE_DIR, "loso_validation_summary.tex")
export_loso_table(loso_summary, loso_table_path)

loso_json_path = joinpath(OUTPUT_ROOT, "loso_summary.json")
loso_json = Dict(
    "n_results" => length(loso_summary.results),
    "mean_validation_rmse" => finite_or_missing(loso_summary.mean_validation_rmse),
    "cross_campaign_stability_score" => finite_or_missing(loso_summary.cross_campaign_stability_score),
    "results" => [
        Dict(
            "val_site_name" => r.val_site_name,
            "discovered_model" => r.discovered_model,
            "coefficients" => r.coefficients,
            "train_rss" => finite_or_missing(r.train_rss),
            "val_rmse" => finite_or_missing(r.val_rmse),
            "val_mae" => finite_or_missing(r.val_mae),
            "val_r2" => finite_or_missing(r.val_r2),
            "coverage_probability" => finite_or_missing(r.coverage_probability),
        ) for r in loso_summary.results
    ],
)
open(loso_json_path, "w") do io
    JSON3.pretty(io, loso_json)
end

summary_data = Dict(
    "generated_at" => string(Dates.now()),
    "n_campaigns" => length(CAMPAIGNS),
    "configuration" => Dict(
        "n_modes" => N_MODES_BENCHMARK,
        "nonlinear_scale" => NONLINEAR_SCALE_BENCHMARK,
        "duration_seconds" => 12.0 * 3600.0,
        "metric_reference" => "final modal state versus campaign-derived observed modal profile",
    ),
    "campaigns" => summary_campaigns,
    "loso" => Dict(
        "n_results" => length(loso_summary.results),
        "mean_validation_rmse" => finite_or_missing(loso_summary.mean_validation_rmse),
        "cross_campaign_stability_score" => finite_or_missing(loso_summary.cross_campaign_stability_score),
    ),
)
open(SUMMARY_PATH, "w") do io
    JSON3.pretty(io, summary_data)
end

p = plot(layout = (1, 2), size = (1200, 500), bottom_margin = 12mm, left_margin = 8mm)
for (j, campaign) in enumerate(sort(collect(keys(profiles))))
    prof = profiles[campaign]
    mode_idx = collect(1:length(prof.observed))
    plot!(
        p[j],
        mode_idx,
        prof.observed,
        lw = 2,
        label = "Observed modal profile",
        xlabel = "Mode Index",
        ylabel = "Amplitude",
        title = "$(campaign) Modal Profile",
    )
    plot!(p[j], mode_idx, prof.final_base, lw = 3.5, ls = :dash, alpha = 0.6, label = "Neutral baseline")
    plot!(p[j], mode_idx, prof.final_phys, lw = 2, label = "Physical closure")
end

fig_path = joinpath(FIG_DIR, "pde_profile_comparison.png")
savefig(p, fig_path)

println("PDE closure benchmark complete.")
println("Summary CSV: ", summary_path)
println("Summary JSON: ", SUMMARY_PATH)
println("Profile figure: ", fig_path)
println("Cross-campaign LaTeX table: ", site_table_path)
println("LOSO LaTeX table: ", loso_table_path)
