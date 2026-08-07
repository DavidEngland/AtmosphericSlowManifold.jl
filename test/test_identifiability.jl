@testset "Identifiability diagnostics" begin
    G = [
        1.0 1.01 0.2;
        2.0 2.02 1.1;
        3.0 3.02 0.1;
        4.0 4.01 2.0;
        5.0 5.02 0.4;
        6.0 6.03 2.8;
    ]

    gram = compute_gram_matrix(G)
    @test size(gram) == (3, 3)
    @test isapprox(gram, transpose(gram); atol = 1e-10)

    mu = compute_mutual_coherence(G)
    @test mu > 0.95

    kappa = compute_condition_number(G)
    @test isfinite(kappa)
    @test kappa > 1.0

    vif = compute_vif(G)
    @test length(vif) == 3
    @test maximum(vif) > 10.0

    cov = compute_parameter_covariance(G; sigma2 = 0.5, ridge = 1e-8)
    @test size(cov) == (3, 3)
    @test isapprox(cov, transpose(cov); atol = 1e-10)

    report = analyze_identifiability(G; sigma2 = 0.5, ridge = 1e-8)
    @test report isa IdentifiabilityReport{Float64}
    @test report.mutual_coherence > 0.95

    pruned = prune_by_mutual_coherence(G; threshold = 0.90)
    @test !isempty(pruned.kept)
    @test !isempty(pruned.dropped)
    @test length(pruned.kept) + length(pruned.dropped) == size(G, 2)
end
