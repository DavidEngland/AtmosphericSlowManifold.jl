function build_pde_system(
    closure::AbstractClosure;
    z_top::Float64 = 3000.0,
    t_end::Float64 = 86400.0,
    coriolis::Float64 = 1e-4,
    v_geostrophic::Float64 = 8.0,
    radiation::Float64 = 0.0
)
    ModelingToolkit.@parameters t z
    ModelingToolkit.@variables u(..) v(..) theta(..)

    Dt = ModelingToolkit.Differential(t)
    Dz = ModelingToolkit.Differential(z)

    state = ManifoldState()
    km = eddy_momentum(closure, state)
    kh = eddy_heat(closure, state)
    flux0 = default_surface_flux(closure, state)

    eqs = [
        Dt(u(t, z)) ~ coriolis * (v(t, z) - v_geostrophic) - Dz(-km * Dz(u(t, z))),
        Dt(v(t, z)) ~ -coriolis * u(t, z) - Dz(-km * Dz(v(t, z))),
        Dt(theta(t, z)) ~ -Dz(-kh * Dz(theta(t, z))) + radiation,
    ]

    bcs = [
        Dz(u(t, 0.0)) ~ flux0,
        Dz(v(t, 0.0)) ~ 0.0,
        Dz(theta(t, 0.0)) ~ 0.0,
        u(t, z_top) ~ v_geostrophic,
        v(t, z_top) ~ 0.0,
        theta(t, z_top) ~ 0.0,
        u(0.0, z) ~ 0.0,
        v(0.0, z) ~ 0.0,
        theta(0.0, z) ~ 0.0,
    ]

    domains = [
        t in Interval(0.0, t_end),
        z in Interval(0.0, z_top),
    ]

    ModelingToolkit.@named pde = ModelingToolkit.PDESystem(eqs, bcs, domains, [t, z], [u(t, z), v(t, z), theta(t, z)])
    return pde
end
