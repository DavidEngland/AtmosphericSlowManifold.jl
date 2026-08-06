@testset "WSINDy closure interface" begin
    ms = ManifoldState()
    km = 0.41 * ms.u_star * ms.z0
    kh = km / 0.74
    flux = ms.u_star^2

    c = WSINDyClosure(km, kh, flux)

    @test isequal(eddy_momentum(c, ms), km)
    @test isequal(eddy_heat(c, ms), kh)
    @test isequal(surface_flux(c, ms), flux)

    v = verify_closure(c)
    @test haskey(v, :km_expr)
    @test haskey(v, :kh_expr)
end
