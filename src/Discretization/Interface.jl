# src/Discretization/Interface.jl
abstract type AbstractDiscretization end

function dispatch_solve(::AbstractDiscretization, pde_sys::PDESystem, closure::AbstractClosure, tspan::Tuple{Float64, Float64}; kwargs...)
    throw(MethodError(dispatch_solve, (AbstractDiscretization, typeof(pde_sys), typeof(closure), typeof(tspan))))
end

function solve_scm(
    pde_sys::PDESystem,
    closure::AbstractClosure,
    disc::AbstractDiscretization,
    tspan::Tuple{Float64, Float64};
    kwargs...
)
    return dispatch_solve(disc, pde_sys, closure, tspan; kwargs...)
end
