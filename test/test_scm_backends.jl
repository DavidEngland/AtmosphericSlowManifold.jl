using LinearAlgebra

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
    @test sg.enable_nonlinear
    @test sg.nonlinear_scale == 1.0
    @test sg.advection_response_scale == 1.0
    @test sg.diffusivity_response_scale == 1.0

    sg_pos = SpectralBLGalerkin(4, 0.75, 100.0)
    @test sg_pos.n_modes == 4
    @test sg_pos.enable_nonlinear
    @test sg_pos.advection_response_scale == 1.0
    @test sg_pos.diffusivity_response_scale == 1.0

    sg_linear = SpectralBLGalerkin(n_modes = 4, lambda = 0.75, H = 100.0, enable_nonlinear = false)
    @test !sg_linear.enable_nonlinear

    sg_adv_only = SpectralBLGalerkin(
        n_modes = 4,
        lambda = 0.75,
        H = 100.0,
        enable_nonlinear = true,
        nonlinear_scale = 1.0,
        advection_response_scale = 1.0,
        diffusivity_response_scale = 0.0,
    )
    sg_diff_only = SpectralBLGalerkin(
        n_modes = 4,
        lambda = 0.75,
        H = 100.0,
        enable_nonlinear = true,
        nonlinear_scale = 1.0,
        advection_response_scale = 0.0,
        diffusivity_response_scale = 1.0,
    )
    @test sg_adv_only.advection_response_scale == 1.0
    @test sg_adv_only.diffusivity_response_scale == 0.0
    @test sg_diff_only.advection_response_scale == 0.0
    @test sg_diff_only.diffusivity_response_scale == 1.0

    tensors = precompute_nonlinear_tensors(sg; n_quad = 96)
    @test tensors isa SpectralNonlinearTensors
    @test size(tensors.triple) == (4, 4, 4)
    @test size(tensors.advection) == (4, 4, 4)
    @test size(tensors.diffusion_flux) == (4, 4, 4)
    @test isfinite(sum(tensors.triple))
    @test isfinite(sum(tensors.advection))
    @test isfinite(sum(tensors.diffusion_flux))
    @test tensors.triple[1, 1, 1] > 0.0

    # Triple-product symmetry in multiplicative indices: C_{kij} = C_{kji}.
    for k in 1:4, i in 1:4, j in 1:4
        @test isapprox(tensors.triple[k, i, j], tensors.triple[k, j, i]; atol = 1e-10)
    end

    a0 = [1.0, -0.3, 0.2, 0.1]
    adv_rhs = AtmosphericSlowManifold._spectral_nonlinear_advection_rhs(tensors, a0, 0.5)
    diff_rhs = AtmosphericSlowManifold._spectral_nonlinear_diffusion_rhs(tensors, a0, 0.7)
    nl_rhs = AtmosphericSlowManifold._spectral_nonlinear_rhs(tensors, a0, 0.5, 0.7)
    @test length(adv_rhs) == 4
    @test length(diff_rhs) == 4
    @test length(nl_rhs) == 4
    @test isapprox(LinearAlgebra.norm(nl_rhs), LinearAlgebra.norm(adv_rhs .+ diff_rhs); atol = 1e-10)
    @test LinearAlgebra.norm(nl_rhs) > 0.0

    adv_only_rhs = AtmosphericSlowManifold._spectral_nonlinear_rhs(tensors, a0, 0.5, 0.0)
    diff_only_rhs = AtmosphericSlowManifold._spectral_nonlinear_rhs(tensors, a0, 0.0, 0.7)
    @test LinearAlgebra.norm(adv_only_rhs) > 0.0
    @test LinearAlgebra.norm(diff_only_rhs) > 0.0
    @test LinearAlgebra.norm(adv_only_rhs .- diff_only_rhs) > 0.0

    # Nonlinear scale-off check.
    zero_rhs = AtmosphericSlowManifold._spectral_nonlinear_rhs(tensors, a0, 0.0, 0.0)
    @test isapprox(LinearAlgebra.norm(zero_rhs), 0.0; atol = 1e-12)

    # Optional solve matrix for local/CI runtime checks.
    if get(ENV, "ASM_RUN_SMOKE", "0") == "1"
        sg_smoke_linear = SpectralBLGalerkin(n_modes = 4, lambda = 0.75, H = 100.0, enable_nonlinear = false)
        sg_smoke_nonlinear = SpectralBLGalerkin(n_modes = 4, lambda = 0.75, H = 100.0, enable_nonlinear = true, nonlinear_scale = 1.0)

        sol_fd_wsindy = solve_scm(pde_wsindy, wsindy, fd, (0.0, 120.0))
        sol_fd_most = solve_scm(pde_most, most, fd, (0.0, 120.0))
        sol_spec_wsindy = solve_scm(pde_wsindy, wsindy, sg_smoke_nonlinear, (0.0, 120.0))
        sol_spec_linear = solve_scm(pde_wsindy, wsindy, sg_smoke_linear, (0.0, 120.0))

        @test string(sol_fd_wsindy.retcode) == "Success"
        @test string(sol_fd_most.retcode) == "Success"
        @test string(sol_spec_wsindy.retcode) == "Success"
        @test string(sol_spec_linear.retcode) == "Success"

        # Nonlinear terms should perturb modal trajectories relative to linear-only evolution.
        final_nl = Vector(sol_spec_wsindy.u[end])
        final_lin = Vector(sol_spec_linear.u[end])
        @test length(final_nl) == length(final_lin)
        @test LinearAlgebra.norm(final_nl .- final_lin) > 1e-8
    end
end

@testset "SpectralBLGalerkin Modal Budget Diagnostics" begin
    disc = SpectralBLGalerkin(6, 0.5, 0.0, 1000.0, true, 1.0, 0.8, 1.2)
    tensors = precompute_nonlinear_tensors(disc)

    u_hat = [1.0, 0.3, -0.1, 0.02, 0.0, 0.0]
    K_hat = [0.15, 0.03, 0.0, 0.0, 0.0, 0.0]

    budget = evaluate_modal_budget(u_hat, K_hat, disc, tensors)

    @test length(budget.linear) == 6
    @test length(budget.advection) == 6
    @test length(budget.diffusion) == 6
    @test length(budget.total) == 6

    expected_total = budget.linear .- (1.0 * 0.8 .* budget.advection) .- (1.0 * 1.2 .* budget.diffusion)
    @test budget.total ≈ expected_total
end
