using Test
using AtmosphericSlowManifold

@testset "Z0HRClosure Unit Tests" begin
    @testset "Constructors and Type Stability" begin
        c = Z0HRClosure()
        @test c.kappa == 0.40
        @test c.pr_t == 1.0
        @test c.B_um == 16.0
        @test c.B_uh == 16.0
        @test c.Ri_c == 0.25
        @test c.eps == 1e-3

        c_f32 = Z0HRClosure(kappa=0.4f0, pr_t=1.0f0, B_um=16.0f0, B_uh=16.0f0, Ri_c=0.25f0, eps=1f-3)
        @test c_f32 isa Z0HRClosure{Float32}
    end

    @testset "Mathematical & Asymptotic Properties" begin
        c = Z0HRClosure()

        # 1. Structural Radicand Positivity: 1 - B_u * Ri_- > 1.0 everywhere
        for Ri in range(-10.0, 10.0, length=100)
            Ri_neg = smooth_min(Ri; eps=c.eps)
            @test (1.0 - c.B_um * Ri_neg) >= 1.0
        end

        # 2. Near-Neutral Point Evaluation (Ri = 0)
        Sm_0, Sh_0 = z0hr_stability_functions(c, 0.0)
        @test isapprox(Sm_0, 0.996; atol=1e-2)
        @test isapprox(Sh_0, 0.996; atol=1e-2)

        # 3. Smooth Physical Cutoff at Critical Richardson Number (Ri >= Ri_c)
        Sm_crit, Sh_crit = z0hr_stability_functions(c, c.Ri_c)
        @test Sm_crit < 1e-4
        @test Sh_crit < 1e-4

        # 4. Asymptotic Convective Regime (Ri -> -infinity)
        _, Sh_unstable = z0hr_stability_functions(c, -5.0)
        @test Sh_unstable > 1.0
    end

    @testset "Manifold State Interfaces" begin
        c = Z0HRClosure()
        state = ManifoldState(z=10.0, z0=0.1, u_star=0.3, r=-0.1)

        Km = eddy_momentum(c, state)
        Kh = eddy_heat(c, state)

        @test Km > 0.0
        @test Kh > 0.0
        @test surface_flux(c, state) == 0.3^2
    end

    @testset "Profile Evaluation Mutators" begin
        c = Z0HRClosure()
        z_grid = collect(range(0.0, 100.0, length=50))
        r_prof = fill(0.02, 50) # Stable stratification

        K_m = zeros(Float64, 50)
        K_h = zeros(Float64, 50)

        evaluate_diffusivity_profile!(K_m, c, z_grid; r_profile=r_prof, u_star=0.3, z0=0.1, mask_below_z0=true)
        evaluate_heat_diffusivity_profile!(K_h, c, z_grid; r_profile=r_prof, u_star=0.3, z0=0.1, mask_below_z0=true)

        @test K_m[1] == 0.0 # Masked at/below z0
        @test all(K_m[2:end] .> 0.0)
        @test all(K_h[2:end] .> 0.0)
    end
end