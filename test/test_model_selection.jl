using Test
using AtmosphericSlowManifold

@testset "Model selection" begin
    m1 = DiscoveredModel{Float64}(:u, OperatorTerm{Float64}[], sqrt(10.0), 1)
    m2 = DiscoveredModel{Float64}(:u, OperatorTerm{Float64}[], sqrt(8.0), 2)
    m3 = DiscoveredModel{Float64}(:u, OperatorTerm{Float64}[], sqrt(12.0), 2)
    m4 = DiscoveredModel{Float64}(:u, OperatorTerm{Float64}[], sqrt(6.0), 3)
    m5 = DiscoveredModel{Float64}(:u, OperatorTerm{Float64}[], sqrt(15.0), 1)
    models = [m1, m2, m3, m4, m5]

    n = 100
    @test isapprox(aic(2, 8.0, n), 2 * 2 + n * log(8.0 / n); atol = 1e-12)
    @test isapprox(bic(2, 8.0, n), 2 * log(n) + n * log(8.0 / n); atol = 1e-12)
    @test isapprox(model_aic(m2, n), aic(m2.sparsity_level, m2.residual_norm^2, n); atol = 1e-12)
    @test isapprox(model_bic(m2, n), bic(m2.sparsity_level, m2.residual_norm^2, n); atol = 1e-12)

    pf_rss = compute_pareto_front(models; objective = :rss)
    @test pf_rss.indices == [1, 2, 4]
    @test pf_rss.complexities == [1, 2, 3]
    @test all(isapprox.(pf_rss.scores, [10.0, 8.0, 6.0]; atol = 1e-10))

    r2_values = [0.60, 0.75, 0.70, 0.80, 0.55]
    pf_r2 = compute_pareto_front(models; objective = :r2, r2_values = r2_values)
    @test pf_r2.indices == [1, 2, 4]
    @test pf_r2.scores == [0.60, 0.75, 0.80]

    @test_throws ArgumentError compute_pareto_front(models; objective = :r2)
    @test_throws ArgumentError compute_pareto_front(models; objective = :bad_objective)

    blocks = ObservationTable[
        ObservationTable(Dict(:u => [1.0, 2.0], :z => [0.0, 1.0]), Dict(:u => "m s^-1", :z => "m")),
        ObservationTable(Dict(:u => [1.0, 2.0, 3.0], :z => [0.0, 1.0, 2.0]), Dict(:u => "m s^-1", :z => "m")),
        ObservationTable(Dict(:u => [1.0, 2.0, 3.0, 4.0], :z => [0.0, 1.0, 2.0, 3.0]), Dict(:u => "m s^-1", :z => "m")),
        ObservationTable(Dict(:u => [1.0, 2.0, 3.0, 4.0, 5.0], :z => [0.0, 1.0, 2.0, 3.0, 4.0]), Dict(:u => "m s^-1", :z => "m")),
        ObservationTable(Dict(:u => [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], :z => [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]), Dict(:u => "m s^-1", :z => "m")),
    ]

    scorer = (model, train_blocks, val_blocks) -> begin
        vlen = sum(length(b.columns[:u]) for b in val_blocks)
        tlen = sum(length(b.columns[:u]) for b in train_blocks)
        return model.residual_norm + 0.01 * tlen + vlen
    end

    cv = kfold_cv_residual(m1, blocks, scorer; k = 2)
    @test cv.k == 2
    @test length(cv.fold_scores) == 2
    @test isfinite(cv.mean_score)

    @test_throws ArgumentError kfold_cv_residual(m1, ObservationTable[], scorer; k = 3)
    @test_throws ArgumentError kfold_cv_residual(m1, blocks, (args...) -> NaN; k = 3)
end
