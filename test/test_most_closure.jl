using Test
using AtmosphericSlowManifold
using Symbolics

function _numeric_manifold_state(; r::Float64, u_star::Float64, z::Float64, z0::Float64 = 0.1)
    return ManifoldState(
        Num(0.0),
        Num(0.0),
        Num(0.0),
        Num(r),
        Num(0.0),
        Num(0.0),
        Num(0.0),
        Num(1.0),
        Num(0.0),
        Num(0.0),
        Num(0.0),
        Num(0.0),
        Num(u_star),
        Num(z),
        Num(z0),
    )
end

@testset "MOST Closure Smooth Blending & Safety Regression Tests" begin
    closure = MOSTClosure()

    # 1. Regression: strong stable stratification (r > 1/gamma) should not throw.
    m_strong_stable = _numeric_manifold_state(r = 0.10, u_star = 0.3, z = 10.0)
    m_extreme_stable = _numeric_manifold_state(r = 1.00, u_star = 0.3, z = 10.0)

    @test_nowarn eddy_momentum(closure, m_strong_stable)
    @test_nowarn eddy_heat(closure, m_strong_stable)
    @test_nowarn eddy_momentum(closure, m_extreme_stable)
    @test_nowarn eddy_heat(closure, m_extreme_stable)

    Km_strong = eddy_momentum(closure, m_strong_stable)
    Kh_strong = eddy_heat(closure, m_strong_stable)
    @test Km_strong > 0.0
    @test Kh_strong > 0.0

    # 2. Smoothness/continuity check near r = 0 transition.
    r_eps = 1e-6
    m_plus = _numeric_manifold_state(r = r_eps, u_star = 0.3, z = 10.0)
    m_minus = _numeric_manifold_state(r = -r_eps, u_star = 0.3, z = 10.0)

    Km_plus = eddy_momentum(closure, m_plus)
    Km_minus = eddy_momentum(closure, m_minus)
    @test isapprox(Km_plus, Km_minus; atol = 1e-4)

    # 3. In-place profile evaluation should stay allocation-free after warmup.
    z_grid = collect(range(0.0, 100.0; length = 64))
    r_prof = fill(0.12, 64)
    K_out = zeros(Float64, 64)

    evaluate_diffusivity_profile!(K_out, closure, z_grid; r_profile = r_prof, z0 = 0.1)

    allocs = @allocated evaluate_diffusivity_profile!(
        K_out,
        closure,
        z_grid;
        r_profile = r_prof,
        z0 = 0.1,
    )
    @test allocs == 0
end