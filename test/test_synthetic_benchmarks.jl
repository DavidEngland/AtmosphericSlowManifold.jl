using Random

@testset "Synthetic benchmarks" begin
    z = collect(range(0.0, 1.0; length = 33))
    t = collect(range(0.0, 2.0; length = 25))

    lin = generate_manufactured_field(z, t; model = :linear, lambda = 0.4)
    nonlin = generate_manufactured_field(z, t; model = :nonlinear, U0 = 1.2, z0 = 0.1, p = 0.25, gamma = 0.2)

    @test size(lin.field) == (length(z), length(t))
    @test size(nonlin.field) == (length(z), length(t))
    @test lin.model == :linear
    @test nonlin.model == :nonlinear

    rng = MersenneTwister(17)
    g = additive_gaussian_noise(lin.field; sigma_rel = 0.05, rng = rng)
    a = ar1_temporal_noise(lin.field; sigma_rel = 0.05, rho = 0.8, rng = rng)
    m = multiplicative_sensor_noise(lin.field; sigma_rel = 0.05, rng = rng)

    @test size(g) == size(lin.field)
    @test size(a) == size(lin.field)
    @test size(m) == size(lin.field)
    @test !all(iszero, g .- lin.field)
    @test !all(iszero, a .- lin.field)
    @test !all(iszero, m .- lin.field)

    miss = apply_missing_data(lin.field; p_drop = 0.2, rng = rng)
    @test size(miss) == size(lin.field)
    @test any(ismissing, miss)

    truth = [1.0, -0.4, 0.0, 2.0]
    direction = [0.2, -0.1, 0.3, -0.4]
    levels = [0.01, 0.05, 0.10, 0.20]
    errs = Float64[]
    for lvl in levels
        est = truth .+ lvl .* direction
        push!(errs, coefficient_l2_error(est, truth))
    end
    @test issorted(errs)

    est_sparse = [1.0, 0.0, 0.1, 1.9]
    pr = precision_recall(est_sparse, truth; tol = 1e-8)
    @test pr.tp == 2
    @test pr.fp == 1
    @test pr.fn == 1
    @test structural_hamming_distance(est_sparse, truth; tol = 1e-8) == 2
    @test false_discovery_rate(est_sparse, truth; tol = 1e-8) > 0.0
    @test equation_sparsity(est_sparse; tol = 1e-8) == 3

    noisy = additive_gaussian_noise(lin.field; sigma_rel = 0.1, rng = MersenneTwister(44))
    @test isfinite(snr_db(lin.field, noisy))

    metrics = evaluate_recovery_metrics(est_sparse, truth)
    @test haskey(Dict(pairs(metrics)), :l2_error)
    @test haskey(Dict(pairs(metrics)), :precision)
    @test haskey(Dict(pairs(metrics)), :recall)
end
