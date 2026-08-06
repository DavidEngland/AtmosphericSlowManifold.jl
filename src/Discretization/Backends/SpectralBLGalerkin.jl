struct SpectralBLGalerkin <: AbstractDiscretization
    n_modes::Int
    lambda::Float64
    H::Float64
end

function SpectralBLGalerkin(; n_modes::Int = 12, lambda::Float64 = 0.75, H::Float64 = 3000.0)
    return SpectralBLGalerkin(n_modes, lambda, H)
end

function dispatch_solve(
    disc::SpectralBLGalerkin,
    pde_sys::PDESystem,
    closure::AbstractClosure,
    tspan::Tuple{Float64, Float64};
    solver = Rodas5P(),
    kwargs...
)
    ModelingToolkit.@variables t
    ModelingToolkit.@variables a(t)[1:disc.n_modes]

    eqs = [ModelingToolkit.Differential(t)(a[n]) ~ -0.1 * n * a[n] for n in 1:disc.n_modes]
    ModelingToolkit.@named modal_ode = ModelingToolkit.ODESystem(eqs, t)
    sys = structural_simplify(modal_ode)

    prob = ODEProblem(sys, ones(disc.n_modes), tspan)
    return solve(prob, solver; kwargs...)
end
