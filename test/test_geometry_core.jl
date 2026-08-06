@testset "Geometry core" begin
    # Fast subsystem: f(x,y,0)=0 => y = x on critical manifold.
    f_fast(x, y, eps; params = Dict{Symbol, Float64}()) = [y[1] - x[1], y[2] - 2x[1]]
    g_slow(x, y, eps; params = Dict{Symbol, Float64}()) = [1.0, -0.5]

    x_path = [[x] for x in range(-1.0, 1.0; length = 5)]
    y_seed = [0.0, 0.0]

    surf = Geometry.solve_critical_surface(f_fast, x_path, y_seed)
    @test length(surf.points) == length(x_path)
    @test maximum(p -> p.residual_norm, surf.points) < 1e-6

    p0 = Geometry.find_manifold_point(f_fast, [0.5], [0.0, 0.0])
    @test isapprox(p0.y[1], 0.5; atol = 1e-6)
    @test isapprox(p0.y[2], 1.0; atol = 1e-6)

    σ = Geometry.fold_indicator(f_fast, p0.x, p0.y)
    @test isfinite(σ)

    slow = Geometry.slow_flow_vector(f_fast, g_slow, p0.x, p0.y)
    desing = Geometry.desingularized_vector_field(f_fast, g_slow, p0.x, p0.y)
    @test length(slow) == 2
    @test length(desing) == 2

    cls = Geometry.classify_folded_singularity([0.0 1.0; -1.0 0.0])
    @test cls in (:folded_focus, :degenerate)

    report = Geometry.fenichel_metrics([2.0 0.0; 0.0 -3.0])
    @test report.min_abs_real > 0
    @test report.is_hyperbolic
end
