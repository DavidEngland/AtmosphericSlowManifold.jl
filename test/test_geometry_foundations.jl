@testset "Geometry foundations" begin
    @variables y1 y2 x1

    f_fast_sym = Num[y1 - x1, y2 - 2x1]
    cache = Geometry.JacobianCache(f_fast_sym, Num[y1, y2], Num[x1])

    Jnum = Geometry.compute_fast_jacobian(cache, Dict(y1 => 0.5, y2 => 1.0, x1 => 0.5))
    @test size(Jnum) == (2, 2)
    @test isapprox(Jnum[1, 1], 1.0; atol = 1e-10)
    @test isapprox(Jnum[2, 2], 1.0; atol = 1e-10)

    adj = Geometry.compute_adjugate([2.0 3.0; 5.0 7.0])
    @test adj == [7.0 -3.0; -5.0 2.0]

    J = [1.0 0.0; 0.0 0.0]
    T = Geometry.evaluate_tangent_space(J)
    @test size(T, 1) == 2

    residual(u, p) = [u[1]^2 - p]
    alg = Geometry.PseudoArclength(ds = 0.05, max_steps = 10)
    branch = Geometry.continue_manifold(residual, [1.0], 1.0, alg)
    @test length(branch.points) >= 2
    @test length(branch.points) == length(branch.parameter_values)
end