function build_pde_system(
    closure::AbstractClosure;
    z_top::Float64 = 3000.0,
    t_end::Float64 = 86400.0,
    coriolis::Float64 = 1e-4,
    v_geostrophic::Float64 = 8.0,
    radiation::Float64 = 0.0,
    z0_value::Float64 = 0.1,
    u_star_value::Float64 = 0.3
)
    ModelingToolkit.@parameters t z
    ModelingToolkit.@variables u(..) v(..) theta(..)

    Dt = ModelingToolkit.Differential(t)
    Dz = ModelingToolkit.Differential(z)

    state = ManifoldState()

    # Bind intrinsic manifold symbols to prognostic/diagnostic fields so
    # closures discovered in manifold coordinates are executable in PDE space.
    manifold_subs = Dict(
        state.eta1 => u(t, z),
        state.eta2 => v(t, z),
        state.eta3 => theta(t, z),
        state.r => sqrt(u(t, z)^2 + v(t, z)^2 + 1e-12),
        state.omega => atan(v(t, z), u(t, z) + 1e-12),
        state.chi => Dz(theta(t, z)),
        state.pi_g => Dz(u(t, z))^2 + Dz(v(t, z))^2,
        state.lambdamin => one(z),
        state.u => u(t, z),
        state.v => v(t, z),
        state.theta => theta(t, z),
        state.q => 0.0,
        state.u_star => u_star_value,
        state.z => z,
        state.z0 => z0_value,
    )

    km = ModelingToolkit.expand_derivatives(Symbolics.substitute(eddy_momentum(closure, state), manifold_subs))
    kh = ModelingToolkit.expand_derivatives(Symbolics.substitute(eddy_heat(closure, state), manifold_subs))
    flux0 = ModelingToolkit.expand_derivatives(Symbolics.substitute(default_surface_flux(closure, state), manifold_subs))

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
