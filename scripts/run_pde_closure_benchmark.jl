using AtmosphericSlowManifold
using CSV
using DataFrames
using Dates
using JSON3
using Plots
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
const RELTOL_BENCHMARK = 1e-4
const ABSTOL_BENCHMARK = 1e-6
const SAVEAT_BENCHMARK = 300.0
const DT_INITIAL_BENCHMARK = 1.0

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

function rmse(a::Vector{Float64}, b::Vector{Float64})
    length(a) == length(b) || throw(DimensionMismatch("Vectors must match for RMSE."))
    return sqrt(mean((a .- b) .^ 2))
end

function r2_score(pred::Vector{Float64}, truth::Vector{Float64})
    length(pred) == length(truth) || throw(DimensionMismatch("Vectors must match for R2."))
    tbar = mean(truth)
    rss = sum((pred .- truth) .^ 2)
    tss = sum((truth .- tbar) .^ 2)
    return tss <= 1e-12 ? 0.0 : max(0.0, 1.0 - rss / tss)
end

function write_text(path::String, content::String)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, content)
    end
    return path
end

function campaign_symbol(campaign::String)
    return Symbol(replace(lowercase(campaign), "-" => "_"))
end

function model_from_profile(
    target::Symbol,
    coeff_u::Float64,
    coeff_du::Float64,
    residual_norm::Float64,
)
    terms = OperatorTerm{Float64}[
        OperatorTerm{Float64}(coeff_u, BasisOperator[BasisOperator(StateVariable(:u), 1.0)]),
        OperatorTerm{Float64}(coeff_du, BasisOperator[BasisOperator(SpatialDerivative(:u, 1), 1.0)]),
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
    p = scatter(
        complexity,
        score,
        xlabel = "Model complexity (k)",
        ylabel = y_label,
        title = title_text,
        markersize = 7,
        label = "Candidates",
    )
    if !isempty(pareto_idx)
        scatter!(
            p,
            complexity[pareto_idx],
            score[pareto_idx],
            markersize = 9,
            markerstrokewidth = 2,
            markercolor = :red,
            label = "Pareto front",
        )
    end
    for i in eachindex(labels)
        annotate!(p, complexity[i], score[i], text(labels[i], 8, :left))
    end
    savefig(p, out_path)
    return out_path
end

function solve_with_fallback(pde_sys, closure, disc, tspan, u0)
    solver_chain = (
        (name = "Rodas5P", alg = Rodas5P()),
        (name = "RadauIIA5", alg = RadauIIA5()),
        (name = "FBDF", alg = FBDF()),
    )

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
                maxiters = 10^7,
                saveat = SAVEAT_BENCHMARK,
                u0 = u0,
            )

            if string(sol.retcode) == "Success"
                return sol, candidate.name
            end
        catch
            # Continue through fallback chain when a solver is incompatible.
        end
    end

    # Return the final attempt result if all candidates fail.
    final = solve_scm(
        pde_sys,
        closure,
        disc,
        tspan;
        solver = Rodas5P(),
        reltol = RELTOL_BENCHMARK,
        abstol = ABSTOL_BENCHMARK,
        dt = DT_INITIAL_BENCHMARK,
        dtmax = 120.0,
        maxiters = 10^7,
        saveat = SAVEAT_BENCHMARK,
        u0 = u0,
    )
    return final, "Rodas5P(final)"
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

    disc = SpectralBLGalerkin(n_modes = N_MODES_BENCHMARK, lambda = 0.75, H = 100.0, enable_nonlinear = true)
    tspan = (0.0, 12.0 * 3600.0)

    mean_wind = campaign_mean_wind(payload.raw)
    mean_wind = isfinite(mean_wind) ? mean_wind : 1.0
    u0 = fill(mean_wind, disc.n_modes)

    pde_phys = build_pde_system(closure; z_top = 100.0, t_end = tspan[2], coriolis = 0.0, v_geostrophic = 0.0, radiation = 0.0)
    pde_base = build_pde_system(baseline; z_top = 100.0, t_end = tspan[2], coriolis = 0.0, v_geostrophic = 0.0, radiation = 0.0)

    sol_phys, solver_phys = solve_with_fallback(pde_phys, closure, disc, tspan, u0)
    sol_base, solver_base = solve_with_fallback(pde_base, baseline, disc, tspan, u0)

    final_phys = Float64.(sol_phys.u[end])
    final_base = Float64.(sol_base.u[end])

    z_grid = collect(range(2.0, 100.0; length = disc.n_modes))
    K_buffer = zeros(Float64, disc.n_modes)
    evaluate_diffusivity_profile!(K_buffer, closure, z_grid)

    r2_phys = r2_score(final_phys, u0)
    r2_base = r2_score(final_base, u0)

    baseline_model = model_from_profile(:u, mean(final_base), 0.0, rmse(final_base, u0))
    physical_model = model_from_profile(:u, mean(final_phys), maximum(K_buffer) / 1000.0, rmse(final_phys, u0))

    candidates = [baseline_model, physical_model]
    closure_labels = ["baseline", "physical"]
    n_obs = length(u0)

    aic_vals = [model_aic(m, n_obs) for m in candidates]
    bic_vals = [model_bic(m, n_obs) for m in candidates]
    rss_vals = [m.residual_norm^2 for m in candidates]
    r2_vals = [r2_base, r2_phys]

    pareto_rss = compute_pareto_front(candidates; objective = :rss)
    pareto_r2 = compute_pareto_front(candidates; objective = :r2, r2_values = r2_vals)

    best_idx = pareto_rss.indices[1]
    best_model = candidates[best_idx]
    best_label = closure_labels[best_idx]

    return (
        campaign = campaign,
        mean_wind = mean_wind,
        rmse_physical = rmse(final_phys, u0),
        rmse_baseline = rmse(final_base, u0),
        r2_physical = r2_phys,
        r2_baseline = r2_base,
        max_km = maximum(K_buffer),
        solver_physical = solver_phys,
        solver_baseline = solver_base,
        phys_retcode = string(sol_phys.retcode),
        baseline_retcode = string(sol_base.retcode),
        final_phys = final_phys,
        final_base = final_base,
        observed = u0,
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
        mean_wind = result.mean_wind,
        rmse_physical = result.rmse_physical,
        rmse_baseline = result.rmse_baseline,
        r2_physical = result.r2_physical,
        r2_baseline = result.r2_baseline,
        max_km = result.max_km,
        phys_retcode = result.phys_retcode,
        baseline_retcode = result.baseline_retcode,
        best_candidate = result.best_label,
        best_aic = minimum(result.aic_vals),
        best_bic = minimum(result.bic_vals),
    ))

    candidate_rows = Any[]
    for i in eachindex(result.candidate_models)
        push!(candidate_rows, Dict(
            "candidate" => result.closure_labels[i],
            "k" => result.candidate_models[i].sparsity_level,
            "residual_norm" => result.candidate_models[i].residual_norm,
            "rss" => result.rss_vals[i],
            "r2" => result.r2_vals[i],
            "aic" => result.aic_vals[i],
            "bic" => result.bic_vals[i],
            "pareto_rss" => in(i, result.pareto_rss.indices),
            "pareto_r2" => in(i, result.pareto_r2.indices),
        ))
    end

    eq_tex = to_latex(result.best_model)
    write_text(joinpath(TABLE_DIR, "$(camp_sym)_best_equation.tex"), string("\\[\n", eq_tex, "\n\\]\n"))
    write_text(
        joinpath(TABLE_DIR, "$(camp_sym)_best_terms.tex"),
        latex_term_table(result.best_model; r2 = maximum(result.r2_vals), residual_norm = result.best_model.residual_norm),
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
    r2_best_by_campaign[camp_sym] = maximum(result.r2_vals)
    residual_best_by_campaign[camp_sym] = result.best_model.residual_norm

    summary_campaigns[campaign] = Dict(
        "n_candidates" => length(result.candidate_models),
        "n_pareto_rss" => length(result.pareto_rss.indices),
        "n_pareto_r2" => length(result.pareto_r2.indices),
        "best_candidate" => result.best_label,
        "best_model" => Dict(
            "terms" => length(result.best_model.terms),
            "residual_norm" => result.best_model.residual_norm,
            "aic" => minimum(result.aic_vals),
            "bic" => minimum(result.bic_vals),
            "latex" => eq_tex,
        ),
        "candidates" => candidate_rows,
    )

    profiles[campaign] = result
end

summary_df = DataFrame(rows)
summary_path = joinpath(TABLE_DIR, "pde_benchmark_summary.csv")
CSV.write(summary_path, summary_df)

site_table_tex = latex_site_summary_table(
    best_models_by_campaign;
    r2_by_site = r2_best_by_campaign,
    residual_by_site = residual_best_by_campaign,
)
site_table_path = joinpath(TABLE_DIR, "cross_campaign_best_models.tex")
write_text(site_table_path, site_table_tex)

summary_data = Dict(
    "generated_at" => string(Dates.now()),
    "n_campaigns" => length(CAMPAIGNS),
    "campaigns" => summary_campaigns,
)
open(SUMMARY_PATH, "w") do io
    JSON3.pretty(io, summary_data)
end

p = plot(layout = (1, 2), size = (1200, 500))
for (j, campaign) in enumerate(sort(collect(keys(profiles))))
    prof = profiles[campaign]
    mode_idx = collect(1:length(prof.observed))
    plot!(
        p[j],
        mode_idx,
        prof.observed,
        lw = 2,
        label = "Observed mean init",
        xlabel = "Mode Index",
        ylabel = "Amplitude",
        title = "$(campaign) Modal Profile",
    )
    plot!(p[j], mode_idx, prof.final_base, lw = 2, ls = :dash, label = "Neutral baseline")
    plot!(p[j], mode_idx, prof.final_phys, lw = 2, label = "Physical closure")
end

fig_path = joinpath(FIG_DIR, "pde_profile_comparison.png")
savefig(p, fig_path)

println("PDE closure benchmark complete.")
println("Summary CSV: ", summary_path)
println("Summary JSON: ", SUMMARY_PATH)
println("Profile figure: ", fig_path)
println("Cross-campaign LaTeX table: ", site_table_path)
