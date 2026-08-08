using CSV
using DataFrames
using JSON3

const ROOT = joinpath(pwd(), "reports", "generated", "pde_benchmark")
const SUMMARY_JSON = joinpath(ROOT, "benchmark_summary.json")
const SUMMARY_CSV = joinpath(ROOT, "tables", "pde_benchmark_summary.csv")
const EXPECTED_CAMPAIGNS = Set(["CASES-99", "SHEBA"])

function check(condition::Bool, message::String)
    println(condition ? "[ok]   " : "[fail] ", message)
    return condition
end

finite_number(value) = value isa Number && isfinite(Float64(value))

function approximately_equal(a, b; rtol::Float64 = 1e-9, atol::Float64 = 1e-10)
    return finite_number(a) && finite_number(b) && isapprox(Float64(a), Float64(b); rtol = rtol, atol = atol)
end

function validate_candidate(candidate, n_observations::Int, campaign::String)
    ok = true
    label = String(candidate["candidate"])
    prefix = "$(campaign) $(label)"

    for metric in ("residual_norm", "rss", "r2", "aic", "bic")
        ok &= check(finite_number(candidate[metric]), "$(prefix) has finite $(metric)")
    end

    if finite_number(candidate["rss"]) && finite_number(candidate["residual_norm"])
        ok &= check(
            approximately_equal(candidate["rss"], Float64(candidate["residual_norm"])^2),
            "$(prefix) satisfies RSS = residual_norm^2",
        )
    end

    if finite_number(candidate["rss"]) && Float64(candidate["rss"]) > 0
        k = Int(candidate["k"])
        rss = Float64(candidate["rss"])
        expected_aic = 2 * k + n_observations * log(rss / n_observations)
        expected_bic = k * log(n_observations) + n_observations * log(rss / n_observations)
        ok &= check(approximately_equal(candidate["aic"], expected_aic), "$(prefix) AIC matches RSS, k, and n")
        ok &= check(approximately_equal(candidate["bic"], expected_bic), "$(prefix) BIC matches RSS, k, and n")
    end

    return ok
end

function main()
    ok = true
    ok &= check(isfile(SUMMARY_JSON), "benchmark summary JSON exists")
    ok &= check(isfile(SUMMARY_CSV), "benchmark summary CSV exists")
    ok || exit(1)

    summary = JSON3.read(read(SUMMARY_JSON, String))
    campaigns = summary["campaigns"]
    present_campaigns = Set(String.(keys(campaigns)))
    ok &= check(present_campaigns == EXPECTED_CAMPAIGNS, "summary contains expected benchmark campaigns")
    ok &= check(haskey(summary, "configuration"), "summary records benchmark configuration")

    for campaign in sort(collect(EXPECTED_CAMPAIGNS))
        haskey(campaigns, campaign) || continue
        result = campaigns[campaign]
        candidates = result["candidates"]
        n_observations = Int(result["n_observations"])

        ok &= check(n_observations > 0, "$(campaign) has a positive observation count")
        ok &= check(length(candidates) == Int(result["n_candidates"]), "$(campaign) candidate count is consistent")
        ok &= check(Int(result["n_pareto_rss"]) > 0, "$(campaign) has a nonempty RSS Pareto front")
        ok &= check(Int(result["n_pareto_r2"]) > 0, "$(campaign) has a nonempty R2 Pareto front")
        ok &= check(!isempty(String(result["observation_source"])), "$(campaign) records its observation source")

        for candidate in candidates
            ok &= validate_candidate(candidate, n_observations, campaign)
        end

        best_label = String(result["best_candidate"])
        best_matches = [candidate for candidate in candidates if String(candidate["candidate"]) == best_label]
        ok &= check(length(best_matches) == 1, "$(campaign) best candidate label is unique")
        isempty(best_matches) && continue

        selected = only(best_matches)
        best_model = result["best_model"]
        ok &= check(Bool(selected["pareto_rss"]), "$(campaign) best candidate belongs to RSS Pareto front")
        for metric in ("residual_norm", "aic", "bic")
            ok &= check(
                approximately_equal(best_model[metric], selected[metric]),
                "$(campaign) best-model $(metric) matches selected candidate",
            )
        end
    end

    table = CSV.read(SUMMARY_CSV, DataFrame)
    required_columns = Set([
        :campaign,
        :observation_source,
        :phys_retcode,
        :baseline_retcode,
        :best_candidate,
        :best_aic,
        :best_bic,
    ])
    ok &= check(required_columns <= Set(Symbol.(names(table))), "benchmark CSV has required result columns")
    ok &= check(Set(String.(table.campaign)) == EXPECTED_CAMPAIGNS, "benchmark CSV contains expected campaigns")
    ok &= check(all(table.phys_retcode .== "Success"), "all physical closure solves succeeded")
    ok &= check(all(table.baseline_retcode .== "Success"), "all baseline solves succeeded")
    ok &= check(all(isfinite, table.best_aic), "all selected AIC values are finite")
    ok &= check(all(isfinite, table.best_bic), "all selected BIC values are finite")

    if ok
        println("PDE benchmark validation passed.")
        return
    end

    println("PDE benchmark validation failed.")
    exit(1)
end

main()