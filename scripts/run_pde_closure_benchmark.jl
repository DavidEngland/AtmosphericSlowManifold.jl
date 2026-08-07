using AtmosphericSlowManifold
using CSV
using DataFrames
using Plots
using Statistics
using DifferentialEquations

const OUTPUT_ROOT = joinpath(pwd(), "reports", "generated", "campaign_exports")
const JSON_DIR = joinpath(OUTPUT_ROOT, "json")
const CSV_DIR = joinpath(OUTPUT_ROOT, "csv")
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

    return (
        campaign = campaign,
        mean_wind = mean_wind,
        rmse_physical = rmse(final_phys, u0),
        rmse_baseline = rmse(final_base, u0),
        max_km = maximum(K_buffer),
        solver_physical = solver_phys,
        solver_baseline = solver_base,
        phys_retcode = string(sol_phys.retcode),
        baseline_retcode = string(sol_base.retcode),
        final_phys = final_phys,
        final_base = final_base,
        observed = u0,
    )
end

rows = NamedTuple[]
profiles = Dict{String, NamedTuple}()
for (campaign, payload) in CAMPAIGNS
    result = run_case(campaign, payload)
    push!(rows, (
        campaign = result.campaign,
        mean_wind = result.mean_wind,
        rmse_physical = result.rmse_physical,
        rmse_baseline = result.rmse_baseline,
        max_km = result.max_km,
        phys_retcode = result.phys_retcode,
        baseline_retcode = result.baseline_retcode,
    ))
    profiles[campaign] = result
end

summary_df = DataFrame(rows)
summary_path = joinpath(TABLE_DIR, "pde_benchmark_summary.csv")
CSV.write(summary_path, summary_df)

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
println("Profile figure: ", fig_path)
