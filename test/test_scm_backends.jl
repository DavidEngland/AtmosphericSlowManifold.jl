@testset "SCM backend scaffolds" begin
    ms = ManifoldState()
    wsindy = WSINDyClosure(0.41 * ms.u_star * ms.z0, 0.41 * ms.u_star * ms.z0 / 0.74, ms.u_star^2)
    most = MOSTClosure()

    pde_wsindy = build_pde_system(wsindy; z_top = 100.0, t_end = 60.0)
    pde_most = build_pde_system(most; z_top = 100.0, t_end = 60.0)

    @test nameof(typeof(pde_wsindy)) == :PDESystem
    @test typeof(pde_wsindy) == typeof(pde_most)
    @test length(pde_wsindy.eqs) == length(pde_most.eqs)

    fd = MethodOfLinesFD(N = 16, H = 100.0, alpha = 2.0, order = 2)
    z = generate_stretched_grid(fd)
    @test length(z) == 16
    @test z[1] ≈ 0.0
    @test z[end] ≈ 100.0

    sg = SpectralBLGalerkin(n_modes = 4, lambda = 0.75, H = 100.0)
    @test sg.n_modes == 4

    # Optional solve matrix for local/CI runtime checks.
    if get(ENV, "ASM_RUN_SMOKE", "0") == "1"
        sol_fd_wsindy = solve_scm(pde_wsindy, wsindy, fd, (0.0, 120.0))
        sol_fd_most = solve_scm(pde_most, most, fd, (0.0, 120.0))
        sol_spec_wsindy = solve_scm(pde_wsindy, wsindy, sg, (0.0, 120.0))

        @test string(sol_fd_wsindy.retcode) == "Success"
        @test string(sol_fd_most.retcode) == "Success"
        @test string(sol_spec_wsindy.retcode) == "Success"
    end
end
