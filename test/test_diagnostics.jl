using Test
using AtmosphericSlowManifold
using Statistics

const D = AtmosphericSlowManifold.Diagnostics

@testset "Error Metrics" begin
    y_true = [1.0, 2.0, 3.0, 4.0]
    y_pred = [1.0, 2.5, 2.0, 5.0]
    y_ref = [1.2, 2.2, 2.8, 4.8]

    @test D.rmse(y_pred, y_true) ≈ sqrt((0.0^2 + 0.5^2 + (-1.0)^2 + 1.0^2) / 4)
    @test D.mae(y_pred, y_true) ≈ (0.0 + 0.5 + 1.0 + 1.0) / 4
    @test D.bias(y_pred, y_true) ≈ (0.0 + 0.5 - 1.0 + 1.0) / 4
    @test D.r2(y_pred, y_true) <= 1.0
    @test isfinite(D.nrmse(y_pred, y_true; norm = :std))
    @test isfinite(D.skill_score(y_pred, y_true, y_ref))
    @test isfinite(D.correlation(y_pred, y_true))
    @test isfinite(D.normalized_bias(y_pred, y_true))
    @test D.relative_l2_error(y_pred, y_true) >= 0.0
    res = D.closure_residual([1.0, 2.0], [1.5, 2.0])
    @test res == [0.5, 0.0]
end

@testset "Energy Budget" begin
    Rn = [100.0, 110.0, 90.0]
    H = [40.0, 45.0, 38.0]
    G = [10.0, 11.0, 9.0]
    LE = [45.0, 50.0, 40.0]

    s = D.surface_energy_budget(Rn, H, G, LE)
    @test s.mean_imbalance ≈ mean(Rn .- H .- G .- LE)
    @test s.rms_imbalance >= 0.0
    @test s.max_imbalance >= 0.0

    Km = [0.1, 0.2, 0.15]
    dudz = [1.0, -2.0, 0.5]
    epsm = D.closure_dissipation(Km, dudz)
    @test epsm == Km .* (dudz .^ 2)

    u = [2.0, 1.5, 1.0]
    v = [0.5, 0.25, 0.0]
    z = [0.0, 10.0, 20.0]
    e = D.energy_residual(u, v, z)
    @test e > 0.0
end

@testset "Manifold Metrics" begin
    vf = [1.0, 0.0, 0.0]
    ns = [0.5, 0.5, 0.0]
    t = D.transversality(vf, ns)
    @test t >= 0.0
    @test t <= 1.0

    state = [1.0, 2.0]
    fold = [1.0 2.0 4.0; 2.0 3.0 5.0]
    fd = D.fold_distance(state, fold)
    @test fd ≈ 0.0

    sm_err = D.slow_manifold_error([1.0, 2.0], x -> x)
    @test sm_err ≈ 0.0

    J = [-2.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -0.2]
    hyp = D.normal_hyperbolicity(J, 1)
    @test hyp >= 1.0
end

@testset "Campaign Diagnostics Aggregation" begin
    u_obs = [2.0, 3.0, 4.0, 5.0]
    u_pred = [2.1, 2.9, 4.2, 4.8]
    theta_obs = [290.0, 291.0, 292.0, 293.0]
    theta_pred = [289.8, 291.2, 292.1, 293.3]

    Rn = [120.0, 118.0, 121.0, 119.0]
    H = [45.0, 44.0, 46.0, 45.0]
    G = [12.0, 12.0, 11.0, 12.0]
    LE = [58.0, 57.0, 59.0, 58.0]

    v_fast = [1.0, 0.2, -0.1]
    n_slow = [0.8, 0.1, 0.0]
    state = [1.0, 2.0, 0.5]
    fold_locus = [
        1.0 1.5 2.0;
        2.0 2.5 3.0;
        0.5 0.6 0.9
    ]

    val_rmse_u = D.rmse(u_pred, u_obs)
    val_rmse_th = D.rmse(theta_pred, theta_obs)
    e_summary = D.surface_energy_budget(Rn, H, G, LE)
    tau_val = D.transversality(v_fast, n_slow)
    dist_fold = D.fold_distance(state, fold_locus)

    summary = D.CampaignDiagnostics(
        "CASES-99_Night3",
        length(u_obs),
        val_rmse_u,
        val_rmse_th,
        0.942,
        0.015,
        e_summary.rms_imbalance,
        4.12,
        tau_val,
        0.03,
        0.18,
        dist_fold,
    )

    @test summary.campaign_name == "CASES-99_Night3"
    @test summary.num_samples == 4
    @test summary.rmse_u == val_rmse_u
    @test summary.rmse_theta == val_rmse_th
    @test summary.energy_imbalance_rms == e_summary.rms_imbalance
    @test summary.transversality_mean == tau_val
    @test summary.fold_distance_min == dist_fold
end
