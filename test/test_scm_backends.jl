@testset "SCM backend scaffolds" begin
    ms = ManifoldState()
    c = WSINDyClosure(0.41 * ms.u_star * ms.z0, 0.41 * ms.u_star * ms.z0 / 0.74, ms.u_star^2)

    pde = build_pde_system(c; z_top = 100.0, t_end = 60.0)
    @test nameof(typeof(pde)) == :PDESystem

    fd = MethodOfLinesFD(N = 16, H = 100.0, alpha = 2.0, order = 2)
    z = generate_stretched_grid(fd)
    @test length(z) == 16
    @test z[1] ≈ 0.0
    @test z[end] ≈ 100.0

    sg = SpectralBLGalerkin(n_modes = 4, lambda = 0.75, H = 100.0)
    @test sg.n_modes == 4

    # Optional smoke run for CI/local envs that opt in.
    if get(ENV, "ASM_RUN_SMOKE", "0") == "1"
        sol_modal = dispatch_solve(sg, pde, c, (0.0, 10.0))
        @test length(sol_modal.u) == 4
    end
end
