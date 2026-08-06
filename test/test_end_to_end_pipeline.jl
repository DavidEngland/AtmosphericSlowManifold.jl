using ModelingToolkit

@testset "End-to-end lifecycle pipeline" begin
    zvals = collect(range(0.0, 100.0; length = 41))
    c0 = 0.15
    c1 = 0.02
    uvals = @. c0 + c1 * zvals

    obs = ObservationTable(
        Dict(
            :z => zvals,
            :u => uvals,
            :v => fill(0.0, length(zvals)),
            :theta => fill(280.0, length(zvals)),
            :q => fill(0.002, length(zvals)),
            :u_star => fill(0.35, length(zvals)),
            :Ri => fill(0.12, length(zvals)),
        ),
        Dict(
            :z => "m",
            :u => "m s^-1",
            :v => "m s^-1",
            :theta => "K",
            :q => "kg kg^-1",
            :u_star => "m s^-1",
            :Ri => "1",
        ),
    )

    ModelingToolkit.@variables z
    candidates = Num[1.0, z]

    discovery = discover_closure(
        obs,
        candidates;
        target = :u,
        basis = GegenbauerBasis(n_spatial = 6, n_temporal = 1, lambda = 0.75),
        lambda = 1e-8,
    )

    @test discovery.discovered_model isa DiscoveredModel{Float64}
    @test discovery.closure isa WSINDyClosure
    @test length(discovery.coefficients) == 2

    pde = build_pde_system(discovery.closure; z_top = 100.0, t_end = 20.0)
    @test nameof(typeof(pde)) == :PDESystem

    disc = SpectralBLGalerkin(n_modes = 4, lambda = 0.75, H = 100.0, enable_nonlinear = true)
    sol = solve_scm(pde, discovery.closure, disc, (0.0, 20.0))
    @test string(sol.retcode) == "Success"
    @test !isempty(sol.u)

    bayes = calibrate(discovery.discovered_model, obs; algorithm = BayesianMCMC(40, 1, 0.8, :NUTS))
    @test bayes isa CalibrationResult{BayesianMCMC}
    @test bayes.diagnostics[:status] == :ok
    @test bayes.diagnostics[:engine] in (:turing, :gaussian_approx)
    @test haskey(bayes.parameters, :posterior_mean)

    vi = calibrate(discovery.discovered_model, obs; algorithm = VariationalInference(50, 200, :ADAM))
    @test vi isa CalibrationResult{VariationalInference}
    @test vi.diagnostics[:status] == :ok
    @test haskey(vi.parameters, :variational_mean)
    @test haskey(vi.parameters, :variational_scale)
end
