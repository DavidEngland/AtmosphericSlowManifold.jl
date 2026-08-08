using Dates
using LinearAlgebra

@testset "Surface forcing ingestion and interpolation" begin
    mktemp() do path, io
        write(io, "campaign,timestamp,station,level,depth_cm,G,Ts,SWC\n")
        write(io, "cases_99,1999-10-01T00:00:00.000,station_1,surface,,10.0,20.0,\n")
        write(io, "cases_99,1999-10-01T00:00:10.000,station_1,surface,,,22.0,\n")
        write(io, "cases_99,1999-10-01T00:00:00.000,station_2,surface,,12.0,21.0,\n")
        write(io, "cases_99,1999-10-01T00:00:00.000,station_1,soil,,,18.0,0.3\n")
        close(io)

        @test_throws ArgumentError load_surface_forcing(path, "cases_99")
        forcing = load_surface_forcing(path, "cases_99"; station="station_1")
        @test forcing.t0 == DateTime(1999, 10, 1)
        @test forcing.t_seconds == [0.0, 10.0]
        @test forcing.Ts ≈ [293.15, 295.15]
        @test forcing.G[1] == 10.0
        @test isnan(forcing.G[2])
        @test interp_forcing(forcing, :Ts, -1.0) ≈ 293.15
        @test interp_forcing(forcing, :Ts, 5.0) ≈ 294.15
        @test interp_forcing(forcing, :Ts, 20.0) ≈ 295.15
        kelvin = load_surface_forcing(
            path,
            "cases_99";
            station="station_1",
            source_temperature_unit=:kelvin,
        )
        @test kelvin.Ts == [20.0, 22.0]
        @test_throws ArgumentError interp_forcing(forcing, :SWC, 0.0)
    end
end

@testset "Spectral surface heat forcing" begin
    mode_count = 4
    phi = AtmosphericSlowManifold._modal_basis_at_z(mode_count, 0.75, 2.0, 100.0)
    forcing = SurfaceForcing(DateTime(2000), [0.0, 10.0], [295.0, 295.0], [NaN, NaN])
    state = [293.0, 0.2, -0.1, 0.05]
    tendency = zeros(mode_count)
    transfer_velocity = 0.01
    AtmosphericSlowManifold._add_surface_heat_forcing!(
        tendency,
        state,
        forcing,
        phi,
        transfer_velocity,
        5.0,
    )

    jacobian = zeros(mode_count, mode_count)
    AtmosphericSlowManifold._add_surface_heat_jacobian!(jacobian, phi, transfer_velocity)
    epsilon = 1e-7
    for column in 1:mode_count
        perturbed_state = copy(state)
        perturbed_state[column] += epsilon
        perturbed_tendency = zeros(mode_count)
        AtmosphericSlowManifold._add_surface_heat_forcing!(
            perturbed_tendency,
            perturbed_state,
            forcing,
            phi,
            transfer_velocity,
            5.0,
        )
        @test isapprox(
            (perturbed_tendency - tendency) / epsilon,
            jacobian[:, column];
            rtol=1e-6,
            atol=1e-8,
        )
    end

    closure = MOSTClosure()
    pde = build_pde_system(closure; z_top=100.0, t_end=10.0)
    disc = SpectralBLGalerkin(n_modes=mode_count, lambda=0.75, H=100.0, enable_nonlinear=false)
    initial_state = [293.0, 0.0, 0.0, 0.0]
    solution = solve_scm(
        pde,
        closure,
        disc,
        (0.0, 10.0);
        u0=initial_state,
        surface_forcing=forcing,
        bulk_transfer_coeff=1e-3,
        reference_wind_speed=1.0,
    )
    @test string(solution.retcode) == "Success"
    @test norm(solution.u[end] - initial_state) > 0.0
    @test_throws ArgumentError solve_scm(
        pde,
        closure,
        disc,
        (0.0, 1.0);
        surface_forcing=forcing,
        bulk_transfer_coeff=-1e-3,
    )
end