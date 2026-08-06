@testset "Discovery modules split" begin
    gf = GegenbauerFamily(0.75, 6)
    val = evaluate_test_function(gf, 2, 50.0, 0.0, 100.0)
    dtt = evaluate_dt_test_function(gf, 2, 1.0, 0.0, 10.0)
    dzz = evaluate_dz2_test_function(gf, 2, 50.0, 0.0, 100.0)
    @test isfinite(val)
    @test isfinite(dtt)
    @test isfinite(dzz)

    lib = build_feature_library([:u], [:Ri], 1)
    obs = ObservationTable(
        Dict(
            :z => collect(range(0.0, 100.0; length = 31)),
            :u => collect(range(1.0, 2.0; length = 31)),
            :Ri => fill(0.2, 31),
            :u_star => fill(0.3, 31),
        ),
        Dict(:z => "m", :u => "m s^-1", :Ri => "1", :u_star => "m s^-1"),
    )
    wf = assemble_weak_system(obs, gf, lib; target = :u, n_temporal = 1)
    @test size(wf.G, 2) == length(lib.features)
    @test length(wf.b) == size(wf.G, 1)

    # Solve a sparse problem with inequality constraint x1 >= 0.
    G = [1.0 0.0; 0.0 1.0; 1.0 1.0]
    b = [1.0, 2.0, 3.0]

    xi_ridge = solve_sparse_regression(G, b, STRidge(lambda = 1e-8, threshold = 1e-10, max_iter = 10))
    @test length(xi_ridge) == 2

    A = [1.0 0.0]
    d = [0.0]
    xi_qp = solve_sparse_regression(G, b, ConstrainedQP(lambda = 1e-8); A_ineq = A, b_ineq = d)
    @test length(xi_qp) == 2
    @test xi_qp[1] >= -1e-8

    # Unified discovery entrypoint with constraints + QP optimizer.
    constraints = AbstractPhysicalConstraint[PositivityConstraint(StateVariable(:u))]
    model_qp = discover(obs, lib, constraints, gf, ConstrainedQP(lambda = 1e-6); target_variable = :K_m, target = :u)
    @test model_qp isa DiscoveredModel{Float64}
    @test model_qp.target_variable == :K_m

    # Unified discovery entrypoint without constraints + STRidge optimizer.
    model_ridge = discover(obs, lib, AbstractPhysicalConstraint[], gf, STRidge(lambda = 1e-8, threshold = 1e-8, max_iter = 10); target_variable = :K_h, target = :u)
    @test model_ridge isa DiscoveredModel{Float64}
    @test model_ridge.target_variable == :K_h

    # STRidge with constraints should currently be rejected.
    @test_throws ArgumentError discover(obs, lib, constraints, gf, STRidge(); target = :u)
end
