using ModelingToolkit

@testset "discover_closure orchestrator" begin
    z = collect(range(0.0, 100.0; length = 81))
    k0 = 0.2
    k1 = 0.015
    y = @. k0 + k1 * z

    obs = ObservationTable(
        Dict(
            :z => z,
            :u => y,
            :u_star => fill(0.35, length(z)),
        ),
        Dict(
            :z => "m",
            :u => "m s^-1",
            :u_star => "m s^-1",
        ),
    )

    ModelingToolkit.@variables z
    candidates = Num[1.0, z]

    result = discover_closure(
        obs,
        candidates;
        target = :u,
        basis = GegenbauerBasis(n_spatial = 6, n_temporal = 1, lambda = 0.75),
        lambda = 1e-8,
    )

    @test result.closure isa WSINDyClosure
    @test length(result.coefficients) == 2
    @test result.discovered_model isa DiscoveredModel{Float64}
    @test result.discovered_model.sparsity_level >= 1
    @test isapprox(result.coefficients[1], k0; atol = 5e-2)
    @test isapprox(result.coefficients[2], k1; atol = 5e-3)
end
