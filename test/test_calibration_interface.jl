@testset "Calibration interface" begin
    terms = OperatorTerm{Float64}[
        OperatorTerm{Float64}(0.2, BasisOperator[BasisOperator(StateVariable(:u), 1.0)]),
        OperatorTerm{Float64}(0.01, BasisOperator[BasisOperator(DiagnosticVariable(:Ri), 1.0)]),
    ]
    model = DiscoveredModel{Float64}(:K_m, terms, 0.0, 2)

    obs = ObservationTable(
        Dict(
            :z => collect(range(0.0, 100.0; length = 21)),
            :u => collect(range(1.0, 2.0; length = 21)),
            :Ri => fill(0.15, 21),
            :u_star => fill(0.3, 21),
            :theta => fill(280.0, 21),
            :v => fill(0.0, 21),
            :q => fill(0.002, 21),
        ),
        Dict(
            :z => "m",
            :u => "m s^-1",
            :Ri => "1",
            :u_star => "m s^-1",
            :theta => "K",
            :v => "m s^-1",
            :q => "kg kg^-1",
        ),
    )

    r_mcmc = calibrate(model, obs; algorithm = BayesianMCMC(200, 2, 0.85, :NUTS))
    @test r_mcmc isa CalibrationResult{BayesianMCMC}
    @test haskey(r_mcmc.parameters, :posterior_mean)
    @test haskey(r_mcmc.parameters, :posterior_std)
    @test haskey(r_mcmc.parameters, :c_1)
    @test length(r_mcmc.parameters[:c_1]) == 400
    @test r_mcmc.diagnostics[:status] == :ok
    @test r_mcmc.diagnostics[:backend] == :BayesianMCMC
    @test r_mcmc.diagnostics[:sampler_type] == :NUTS
    @test r_mcmc.diagnostics[:engine] in (:turing, :gaussian_approx)
    @test length(r_mcmc.parameters[:posterior_std]) == length(r_mcmc.parameters[:posterior_mean])

    r_mle = calibrate(model, obs; algorithm = MaximumLikelihood(:LBFGS, 100, 1e-6))
    @test r_mle isa CalibrationResult{MaximumLikelihood}
    @test haskey(r_mle.parameters, :mle)
    @test r_mle.diagnostics[:backend] == :MaximumLikelihood

    r_vi = calibrate(model, obs; algorithm = VariationalInference(100, 300, :ADAM))
    @test r_vi isa CalibrationResult{VariationalInference}
    @test haskey(r_vi.parameters, :variational_mean)
    @test haskey(r_vi.parameters, :variational_scale)
    @test haskey(r_vi.parameters, :c_1)
    @test length(r_vi.parameters[:c_1]) == 100
    @test r_vi.diagnostics[:backend] == :VariationalInference
    @test r_vi.diagnostics[:status] == :ok
    @test haskey(r_vi.diagnostics, :elbo_trace)
    @test !isempty(r_vi.diagnostics[:elbo_trace])

    # Submodule API should mirror root exports.
    r_sub = Calibration.calibrate(model, obs; algorithm = Calibration.BayesianMCMC(1000, 4, 0.8, :NUTS))
    @test r_sub isa Calibration.CalibrationResult{Calibration.BayesianMCMC}
    @test haskey(r_sub.parameters, :c_2)

    # Deterministic recovery benchmark: u = c1*v + c2*Ri.
    vcol = collect(range(-1.0, 1.0; length = 31))
    ricol = collect(range(0.1, 0.4; length = 31))
    c1_true = 2.0
    c2_true = -1.0
    ucol = @. c1_true * vcol + c2_true * ricol

    obs_fit = ObservationTable(
        Dict(
            :z => collect(range(0.0, 100.0; length = 31)),
            :u => ucol,
            :v => vcol,
            :Ri => ricol,
            :u_star => fill(0.3, 31),
            :theta => fill(280.0, 31),
            :q => fill(0.002, 31),
        ),
        Dict(
            :z => "m",
            :u => "m s^-1",
            :v => "m s^-1",
            :Ri => "1",
            :u_star => "m s^-1",
            :theta => "K",
            :q => "kg kg^-1",
        ),
    )

    fit_terms = OperatorTerm{Float64}[
        OperatorTerm{Float64}(0.0, BasisOperator[BasisOperator(StateVariable(:v), 1.0)]),
        OperatorTerm{Float64}(0.0, BasisOperator[BasisOperator(DiagnosticVariable(:Ri), 1.0)]),
    ]
    fit_model = DiscoveredModel{Float64}(:u, fit_terms, 0.0, 2)

    fit = calibrate(fit_model, obs_fit; algorithm = MaximumLikelihood(:LBFGS, 200, 1e-10))
    @test fit isa CalibrationResult{MaximumLikelihood}
    @test isapprox(fit.parameters[:mle][1], c1_true; atol = 1e-5)
    @test isapprox(fit.parameters[:mle][2], c2_true; atol = 1e-5)
    @test fit.diagnostics[:r2] > 0.999999
end
