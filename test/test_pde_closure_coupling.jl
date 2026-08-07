using Test
using AtmosphericSlowManifold
using LinearAlgebra

@testset "PhysicalClosure Struct & Ingestion" begin
    mktempdir() do d
        json_path = joinpath(d, "campaign_model_and_diagnostics.json")
        open(json_path, "w") do io
            write(io, """
            {
              "diagnostics": {
                "obukhov_scaling": {"mean": -80.0, "n": 10, "present": true},
                "similarity_parameters": {
                  "zeta": {"mean": 0.2, "n": 10},
                  "phi_obs": {"mean": 1.3, "n": 10}
                },
                "stats": {
                  "ustar": {"mean": 0.25, "n": 10}
                }
              },
              "terms": [
                {"name": "zeta", "coefficient": 1.0},
                {"name": "phi_obs", "coefficient": 0.5}
              ]
            }
            """)
        end

        c = PhysicalSimilarityClosure(json_path)
        @test c isa PhysicalSimilarityClosure
        @test isfinite(c.ustar)
        @test isfinite(c.L_obukhov)
        @test c.karman > 0.0

        z_grid = collect(range(2.0, 100.0; length = 18))
        K = zeros(Float64, length(z_grid))
        evaluate_diffusivity_profile!(K, c, z_grid)
        alloc = @allocated evaluate_diffusivity_profile!(K, c, z_grid)
        @test alloc == 0
        @test all(isfinite, K)
    end
end

@testset "PhysicalClosure diagnostics filename fallback" begin
    mktempdir() do d
        json_model = joinpath(d, "sheba_model_and_diagnostics.json")
        open(json_model, "w") do io
            write(io, """
            {
              "diagnostics": {
                "obukhov_scaling": {"mean": -50.0, "n": 12, "present": true},
                "similarity_parameters": {
                  "zeta": {"mean": 0.2, "n": 12},
                  "phi_obs": {"mean": 1.1, "n": 12}
                },
                "stats": {
                  "ustar": {"mean": 0.3, "n": 12}
                }
              },
              "terms": [
                {"name": "zeta", "coefficient": 1.0},
                {"name": "phi_obs", "coefficient": 0.25}
              ]
            }
            """)
        end

        # Request *_diagnostics.json even though only *_model_and_diagnostics.json exists.
        c = PhysicalSimilarityClosure(joinpath(d, "sheba_diagnostics.json"))
        @test c isa PhysicalSimilarityClosure
        @test isfinite(c.L_obukhov)
        @test isfinite(c.ustar)
    end
end

@testset "PDE RHS Allocation & Stability" begin
    closure = PhysicalSimilarityClosure(
        phi_coeffs = [1.0, 0.15],
        zeta_coeffs = [0.0, 1.0],
        karman = 0.4,
        ustar = 0.3,
        L_obukhov = -75.0,
        z_ref = 10.0,
    )

    disc = SpectralBLGalerkin(n_modes = 8, lambda = 0.75, H = 100.0, enable_nonlinear = true)
    ws = AtmosphericSlowManifold.build_boundary_layer_workspace(disc, closure)
    AtmosphericSlowManifold.update_diffusivity_buffers!(ws)

    L = AtmosphericSlowManifold._linear_modal_operator(disc)
    tensors = precompute_nonlinear_tensors(disc; n_quad = 96)

    u = ones(Float64, disc.n_modes)
    du = similar(u)

    AtmosphericSlowManifold.spectral_rhs!(du, u, L, tensors, ws, disc, 0.2, 0.1)
    function _run_rhs_loop!(du_loc, u_loc, L_loc, tensors_loc, ws_loc, disc_loc)
      for _ in 1:100
        AtmosphericSlowManifold.spectral_rhs!(du_loc, u_loc, L_loc, tensors_loc, ws_loc, disc_loc, 0.2, 0.1)
        end
    end
    _run_rhs_loop!(du, u, L, tensors, ws, disc)
    alloc = @allocated _run_rhs_loop!(du, u, L, tensors, ws, disc)

    @test alloc == 0
    @test all(isfinite, du)
end

@testset "Coupled Time Integration" begin
    closure = PhysicalSimilarityClosure(
        phi_coeffs = [1.0],
        zeta_coeffs = [0.0, 1.0],
        karman = 0.4,
        ustar = 0.3,
        L_obukhov = -100.0,
        z_ref = 10.0,
    )

    pde = build_pde_system(closure; z_top = 100.0, t_end = 3600.0, coriolis = 0.0, v_geostrophic = 0.0, radiation = 0.0)
    disc = SpectralBLGalerkin(n_modes = 8, lambda = 0.75, H = 100.0, enable_nonlinear = false)

    u0 = collect(range(2.0, 0.5; length = disc.n_modes))
    sol = solve_scm(pde, closure, disc, (0.0, 3600.0); u0 = u0)

    e0 = 0.5 * sum(abs2, Float64.(u0))
    e1 = 0.5 * sum(abs2, Float64.(sol.u[end]))

    @test string(sol.retcode) == "Success"
    @test e1 <= e0 + 1e-6
end
