using CSV
using DataFrames
using JSON3
using Statistics

@testset "LOSO Cross-Validation" begin
    function mk_campaign(name::String, shift::Float64; n::Int = 24)
        x = collect(range(0.0, 1.0; length = n))
        y = 1.0 .+ 2.5 .* x .+ shift
        states = hcat(x, y)
        t = collect(range(0.0, n - 1; length = n))
        z = collect(range(2.0, 100.0; length = n))
        return CampaignData(name, states, t, z)
    end

    # 1) create_loso_splits on 4 synthetic campaigns.
    campaigns4 = CampaignData[
        mk_campaign("SHEBA", 0.00),
        mk_campaign("CASES-99", 0.02),
        mk_campaign("FLOSS", -0.01),
        mk_campaign("BLLAST", 0.01),
    ]

    splits = create_loso_splits(campaigns4)
    @test length(splits) == 4
    @test all(length(s.train_campaigns) == 3 for s in splits)
    @test all(s.val_campaign.site_name == s.val_site_name for s in splits)
    @test sort([s.val_site_name for s in splits]) == ["BLLAST", "CASES-99", "FLOSS", "SHEBA"]

    # 2) evaluate_out_of_sample deterministic metrics.
    X = reshape([0.0, 1.0, 2.0, 3.0], :, 1)
    y_true = [0.0, 1.0, 2.0, 3.0]
    coeff = [0.5]
    metrics = evaluate_out_of_sample(coeff, X, y_true)

    y_pred = X * coeff
    resid = y_pred .- y_true
    expected_rmse = sqrt(mean(abs2, resid))
    expected_mae = mean(abs.(resid))
    expected_r2 = 1.0 - sum(abs2, resid) / sum(abs2, y_true .- mean(y_true))
    sigma = std(resid)
    expected_coverage = mean(abs.(resid) .<= 1.96 * sigma)

    @test isapprox(metrics.rmse, expected_rmse; atol = 1e-12)
    @test isapprox(metrics.mae, expected_mae; atol = 1e-12)
    @test isapprox(metrics.r2, expected_r2; atol = 1e-12)
    @test isapprox(metrics.coverage_probability, expected_coverage; atol = 1e-12)

    # 3) run_loso_cross_validation end-to-end on 3 campaigns.
    campaigns3 = CampaignData[
        mk_campaign("SHEBA", 0.00),
        mk_campaign("CASES-99", 0.03),
        mk_campaign("FLOSS", -0.02),
    ]

    mock_discovery = function (_train_campaigns, Xtrain, ytrain; kwargs...)
        coeffs = Xtrain \ ytrain
        rss = sum(abs2, Xtrain * coeffs .- ytrain)
        return (coefficients = coeffs, discovered_model = "mock_linear", train_rss = rss)
    end

    summary = run_loso_cross_validation(campaigns3, mock_discovery; feature_cols = [1], target_col = 2)

    @test summary isa LOSOSummary
    @test length(summary.results) == 3
    @test summary.mean_validation_rmse >= 0.0
    @test isfinite(summary.cross_campaign_stability_score)

    expected_sites = Set(["SHEBA", "CASES-99", "FLOSS"])
    @test Set([r.val_site_name for r in summary.results]) == expected_sites

    for r in summary.results
        @test r.discovered_model == "mock_linear"
        @test length(r.coefficients) == 1
        @test isfinite(r.train_rss)
        @test r.val_rmse >= 0.0
        @test r.val_mae >= 0.0
        @test isfinite(r.val_r2)
        @test 0.0 <= r.coverage_probability <= 1.0
    end
end

@testset "Campaign Artifact Ingestion" begin
    mktempdir() do tmpdir
        sites = ["SHEBA", "CASES-99", "FLOSS"]

        for site in sites
            csv_file = joinpath(tmpdir, "$(site).csv")
            open(csv_file, "w") do io
                println(io, "time,z,feature1,feature2,target")
                for t in 1:10
                    println(io, "$(t),$(t * 0.5),$(1.0 + t),$(2.0 * t),$(0.5 * t + 1.0)")
                end
            end

            json_file = joinpath(tmpdir, "$(site).json")
            open(json_file, "w") do io
                println(io, "{\"latitude\": 70.0, \"instrument\": \"radiometer\"}")
            end
        end

        sheba = load_campaign_data(tmpdir, "SHEBA")
        @test sheba.site_name == "SHEBA"
        @test size(sheba.states) == (10, 3)
        @test length(sheba.time) == 10
        @test length(sheba.z) == 10
        @test sheba.auxiliary[:latitude] == 70.0

        campaigns = load_all_campaigns(tmpdir; sites = sites)
        @test length(campaigns) == 3
        @test [c.site_name for c in campaigns] == sites

        mock_discovery(c_train, X, y; kwargs...) = (coefficients = X \ y, model = "linear_fit")
        summary = run_artifact_loso(tmpdir, mock_discovery; sites = sites)
        @test summary isa LOSOSummary
        @test length(summary.results) == 3
        @test summary.mean_validation_rmse >= 0.0
        @test isfinite(summary.cross_campaign_stability_score)
    end
end

@testset "LOSO Artifact Ingestion" begin
    @testset "Mock Artifact Roundtrip" begin
        mktempdir() do d
            mkpath(joinpath(d, "csv"))
            mkpath(joinpath(d, "json"))

            U = [1.0 2.0 3.0; 4.0 5.0 6.0]
            z = [10.0, 20.0]
            t = [0.0, 1.0, 2.0]

            CSV.write(joinpath(d, "csv", "sheba_profile.csv"), DataFrame(U, :auto))
            CSV.write(joinpath(d, "csv", "sheba_z.csv"), DataFrame(z = z))
            CSV.write(joinpath(d, "csv", "sheba_time.csv"), DataFrame(t = t))

            open(joinpath(d, "json", "sheba_model_and_diagnostics.json"), "w") do io
                JSON3.pretty(io, Dict("site" => "SHEBA", "n_profiles" => 1, "params" => Dict("alpha" => 2.0)))
            end

            c = load_campaign_data(d, "SHEBA")
            @test c.site_name == "SHEBA"
            @test size(c.states) == (3, 2)
            @test c.time == t
            @test c.z == z
            @test c.auxiliary[:ingestion_mode] == :profile_bundle
            @test haskey(c.auxiliary, :metadata)
            @test haskey(c.auxiliary, :metadata_source)
        end
    end

    @testset "Benchmark Artifact Verification" begin
        repo_root = dirname(@__DIR__)
        artifacts_dir = joinpath(repo_root, "reports", "generated", "campaign_exports")

        if isdir(artifacts_dir)
            campaigns = load_all_campaigns(artifacts_dir; sites = ["SHEBA", "CASES-99", "FLOSS", "BLLAST"])
            @test length(campaigns) >= 2
            for c in campaigns
                @test !isempty(c.time)
                @test !isempty(c.z)
                @test size(c.states, 1) > 0
                @test size(c.states, 2) > 1
            end
        else
            # Keep suite robust in clean environments with no generated artifacts.
            @test true
        end
    end

    @testset "End-to-End Ingestion LOSO" begin
        mktempdir() do d
            mkpath(joinpath(d, "csv"))
            mkpath(joinpath(d, "json"))

            function write_site(site::String, shift::Float64)
                slug = replace(lowercase(site), "-" => "_")
                n = 20
                time = collect(0.0:(n - 1))
                zeta = collect(range(0.0, 1.0; length = n))
                phi = 1.0 .+ 1.5 .* zeta .+ shift
                ustar = 0.2 .+ 0.05 .* zeta
                df = DataFrame(time = time, zeta = zeta, phi_obs = phi, ustar = ustar)
                CSV.write(joinpath(d, "csv", "$(slug)_raw.csv"), df)

                open(joinpath(d, "json", "$(slug)_model_and_diagnostics.json"), "w") do io
                    JSON3.pretty(io, Dict("site" => site, "shift" => shift))
                end
            end

            write_site("SHEBA", 0.00)
            write_site("CASES-99", 0.03)
            write_site("FLOSS", -0.02)

            mock_discovery = function (_train_campaigns, Xtrain, ytrain; kwargs...)
                coeffs = Xtrain \ ytrain
                rss = sum(abs2, Xtrain * coeffs .- ytrain)
                return (coefficients = coeffs, discovered_model = "artifact_mock_linear", train_rss = rss)
            end

            summary = run_artifact_loso(
                d,
                mock_discovery;
                sites = ["SHEBA", "CASES-99", "FLOSS"],
                feature_cols = [1],
                target_col = 2,
            )

            @test summary isa LOSOSummary
            @test length(summary.results) == 3
            @test summary.mean_validation_rmse >= 0.0
            @test isfinite(summary.cross_campaign_stability_score)
        end
    end
end
