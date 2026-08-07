using ModelingToolkit

@testset "WSINDy engine" begin
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
    basis = GegenbauerBasis(n_spatial = 6, n_temporal = 1, lambda = 0.75)

    G, b = build_weak_library(obs, basis, candidates; target = :u)
    coeffs = fit_wsindy_jump(G, b, candidates; lambda = 1e-8, positivity_constraints = true)

    @test length(coeffs) == 2
    @test isapprox(coeffs[1], k0; atol = 5e-2)
    @test isapprox(coeffs[2], k1; atol = 5e-3)

    closure = extract_closure(coeffs, candidates)
    @test closure isa WSINDyClosure

    ms = ManifoldState()
    km_expr = eddy_momentum(closure, ms)
    @test km_expr isa Num
end

@testset "Constrained STLSQ" begin
    G = [
        1.0 0.0 1.0
        0.0 1.0 1.0
        1.0 1.0 2.0
        2.0 1.0 3.0
    ]
    b = [0.5, 0.5, 1.0, 1.5]

    lb = [0.0, 0.0, 0.0]
    beta = constrained_stlsq(G, b, 1e-6; lower_bounds = lb, max_iter = 50, tol = 1e-10)

    @test length(beta) == 3
    @test all(beta .>= -1e-10)

    beta_sparse = constrained_stlsq(G, b, 0.2; lower_bounds = lb, max_iter = 50)
    @test count(abs.(beta_sparse) .> 1e-8) <= count(abs.(beta) .> 1e-8)
end
