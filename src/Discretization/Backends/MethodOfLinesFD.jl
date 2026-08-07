# src/Discretization/Backends/MethodOfLinesFD.jl
struct MethodOfLinesFD <: AbstractDiscretization
    N::Int
    H::Float64
    alpha::Float64
    order::Int
end

function MethodOfLinesFD(; N::Int = 100, H::Float64 = 3000.0, alpha::Float64 = 3.5, order::Int = 2)
    return MethodOfLinesFD(N, H, alpha, order)
end

generate_stretched_grid(d::MethodOfLinesFD) = generate_stretched_grid(d.N, d.H, d.alpha)

function dispatch_solve(
    disc::MethodOfLinesFD,
    pde_sys::PDESystem,
    closure::AbstractClosure,
    tspan::Tuple{Float64, Float64};
    solver = TRBDF2(),
    kwargs...
)
    z_grid = generate_stretched_grid(disc)
    ivs = pde_sys.ivs
    tvar = ivs[1]
    zvar = ivs[2]

    mol_disc = MOLFiniteDifference([zvar => z_grid], tvar; approx_order = disc.order)
    prob = discretize(pde_sys, mol_disc)

    return solve(prob, solver; kwargs...)
end
