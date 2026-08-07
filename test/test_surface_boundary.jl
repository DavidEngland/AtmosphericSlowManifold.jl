using Test
using AtmosphericSlowManifold
using Symbolics

function _surface_test_state(; u_star::Float64 = 0.3, theta::Float64 = 290.0, z::Float64 = 10.0, z0::Float64 = 0.1)
    return ManifoldState(
        Num(0.0),
        Num(0.0),
        Num(0.0),
        Num(0.0),
        Num(0.0),
        Num(0.0),
        Num(0.0),
        Num(1.0),
        Num(0.0),
        Num(0.0),
        Num(theta),
        Num(0.0),
        Num(u_star),
        Num(z),
        Num(z0),
    )
end

@testset "SurfaceBoundary helpers" begin
    c = MOSTClosure()
    m = _surface_test_state()

    tau = default_surface_flux(c, m)
    @test tau == surface_flux(c, m)

    # Current ManifoldState has no q_theta/theta_s, so fallback should be zero.
    qh = default_surface_heat_flux(c, m)
    @test qh == zero(m.u_star)

    bc = surface_boundary_conditions(c, m)
    @test haskey(bc, :momentum_flux)
    @test haskey(bc, :heat_flux)
    @test bc.momentum_flux == tau
    @test bc.heat_flux == qh
end
