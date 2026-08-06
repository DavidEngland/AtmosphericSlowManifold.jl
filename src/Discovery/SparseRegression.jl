using HiGHS

abstract type AbstractSparseOptimizer end

struct STRidge <: AbstractSparseOptimizer
    lambda::Float64
    threshold::Float64
    max_iter::Int
end

struct ConstrainedQP <: AbstractSparseOptimizer
    lambda::Float64
    solver::Any
end

STRidge(; lambda::Float64 = 1e-3, threshold::Float64 = 1e-2, max_iter::Int = 50) =
    STRidge(lambda, threshold, max_iter)

ConstrainedQP(; lambda::Float64 = 1e-3, solver = HiGHS.Optimizer) =
    ConstrainedQP(lambda, solver)

function _solve_ridge(G::Matrix{Float64}, b::Vector{Float64}, lambda::Float64)
    n = size(G, 2)
    A = G' * G + lambda * I(n)
    rhs = G' * b
    return A \ rhs
end

function solve_sparse_regression(G::Matrix{Float64}, b::Vector{Float64}, opt::STRidge)
    size(G, 1) == length(b) || throw(ArgumentError("G and b dimensions are inconsistent."))

    n = size(G, 2)
    support = trues(n)
    xi = zeros(n)

    for _ in 1:opt.max_iter
        if count(support) == 0
            break
        end

        Gs = G[:, support]
        xis = _solve_ridge(Gs, b, opt.lambda)
        xi_new = zeros(n)
        xi_new[support] .= xis

        new_support = abs.(xi_new) .>= opt.threshold
        if new_support == support
            xi = xi_new
            break
        end

        xi = xi_new
        support = new_support
    end

    return xi
end

function solve_sparse_regression(
    G::Matrix{Float64},
    b::Vector{Float64},
    opt::ConstrainedQP;
    A_ineq::Union{Nothing, Matrix{Float64}} = nothing,
    b_ineq::Union{Nothing, Vector{Float64}} = nothing,
)
    size(G, 1) == length(b) || throw(ArgumentError("G and b dimensions are inconsistent."))

    n_features = size(G, 2)
    model = JuMP.Model(opt.solver)
    JuMP.set_silent(model)

    JuMP.@variable(model, xi[1:n_features])
    JuMP.@variable(model, abs_xi[1:n_features] >= 0)

    JuMP.@constraint(model, [i in 1:n_features], xi[i] <= abs_xi[i])
    JuMP.@constraint(model, [i in 1:n_features], -xi[i] <= abs_xi[i])

    if !(A_ineq === nothing)
        b_ineq === nothing && throw(ArgumentError("b_ineq must be provided when A_ineq is specified."))
        size(A_ineq, 2) == n_features || throw(ArgumentError("A_ineq must have $(n_features) columns."))
        size(A_ineq, 1) == length(b_ineq) || throw(ArgumentError("A_ineq row count must match b_ineq length."))

        for i in 1:size(A_ineq, 1)
            JuMP.@constraint(model, sum(A_ineq[i, j] * xi[j] for j in 1:n_features) >= b_ineq[i])
        end
    end

    JuMP.@objective(
        model,
        Min,
        0.5 * sum((sum(G[i, j] * xi[j] for j in 1:n_features) - b[i])^2 for i in 1:size(G, 1)) +
        opt.lambda * sum(abs_xi),
    )

    JuMP.optimize!(model)
    status = JuMP.termination_status(model)
    status in (JuMP.MOI.OPTIMAL, JuMP.MOI.LOCALLY_SOLVED) ||
        throw(ArgumentError("Sparse QP optimization failed with status $(status)."))

    return Float64.(JuMP.value.(xi))
end
