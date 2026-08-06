module Continuation

using LinearAlgebra

export AbstractContinuationAlgorithm, PseudoArclength, ContinuationBranch, continue_manifold, continue_set

abstract type AbstractContinuationAlgorithm end

struct PseudoArclength <: AbstractContinuationAlgorithm
    ds::Float64
    ds_min::Float64
    ds_max::Float64
    max_steps::Int
end

PseudoArclength(; ds::Float64 = 1e-2, ds_min::Float64 = 1e-5, ds_max::Float64 = 1e-1, max_steps::Int = 500) =
    PseudoArclength(ds, ds_min, ds_max, max_steps)

struct ContinuationBranch
    points::Vector{Vector{Float64}}
    tangents::Vector{Vector{Float64}}
    parameter_values::Vector{Float64}
    stability::Vector{Bool}
end

function _fd_jacobian_u(residual_fn, u::Vector{Float64}, p::Float64; eps::Float64 = 1e-6)
    r0 = residual_fn(u, p)
    m = length(r0)
    n = length(u)
    J = zeros(m, n)
    for j in 1:n
        up = copy(u)
        um = copy(u)
        up[j] += eps
        um[j] -= eps
        rp = residual_fn(up, p)
        rm = residual_fn(um, p)
        J[:, j] = (rp .- rm) ./ (2eps)
    end
    return J
end

function _newton_correct(residual_fn, u_init::Vector{Float64}, p::Float64; tol::Float64 = 1e-8, maxiter::Int = 20)
    u = copy(u_init)
    for _ in 1:maxiter
        r = residual_fn(u, p)
        norm(r) <= tol && return u, true
        J = _fd_jacobian_u(residual_fn, u, p)
        size(J, 1) == size(J, 2) || return u, false
        du = J \ r
        u .-= du
        norm(du) <= tol && return u, true
    end
    return u, false
end

function continue_manifold(
    residual_fn,
    u0::Vector{Float64},
    p0::Float64,
    alg::PseudoArclength;
    tangent0::Union{Nothing, Vector{Float64}} = nothing,
    stability_fn::Union{Nothing, Function} = nothing,
)
    points = Vector{Vector{Float64}}()
    tangents = Vector{Vector{Float64}}()
    parameter_values = Float64[]
    stability = Bool[]

    u = copy(u0)
    p = p0
    t = tangent0 === nothing ? ones(length(u0)) : copy(tangent0)
    norm(t) > 0 || (t .= 1.0)
    t ./= norm(t)

    push!(points, copy(u))
    push!(tangents, copy(t))
    push!(parameter_values, p)
    push!(stability, stability_fn === nothing ? true : Bool(stability_fn(u, p)))

    ds = alg.ds
    for _ in 1:alg.max_steps
        u_pred = u .+ ds .* t
        p = p + ds

        u_next, ok = _newton_correct(residual_fn, u_pred, p)
        if !ok
            ds = max(0.5 * ds, alg.ds_min)
            ds <= alg.ds_min && break
            continue
        end

        t_next = u_next .- u
        if norm(t_next) > 0
            t = t_next / norm(t_next)
        end

        u = u_next
        push!(points, copy(u))
        push!(tangents, copy(t))
        push!(parameter_values, p)
        push!(stability, stability_fn === nothing ? true : Bool(stability_fn(u, p)))

        ds = min(1.05 * ds, alg.ds_max)
    end

    return ContinuationBranch(points, tangents, parameter_values, stability)
end

function continue_set(set, param_idx::Int, alg::AbstractContinuationAlgorithm)
    if hasproperty(set, :points) && set.points isa AbstractMatrix
        ncols = size(set.points, 2)
        ncols > 0 || throw(ArgumentError("Invariant set has no points to continue."))
        param_idx <= size(set.points, 1) || throw(ArgumentError("param_idx out of bounds for set point coordinates."))

        p0 = set.points[param_idx, 1]
        u0 = Float64.(set.points[1:(param_idx - 1), 1])
        residual_fn = (u, p) -> u .- u0
        return continue_manifold(residual_fn, u0, p0, alg)
    elseif hasproperty(set, :coordinates) && set.coordinates isa AbstractMatrix
        ncols = size(set.coordinates, 2)
        ncols > 0 || throw(ArgumentError("Invariant set has no coordinates to continue."))
        param_idx <= size(set.coordinates, 1) || throw(ArgumentError("param_idx out of bounds for set coordinates."))

        p0 = set.coordinates[param_idx, 1]
        u0 = Float64.(set.coordinates[1:(param_idx - 1), 1])
        residual_fn = (u, p) -> u .- u0
        return continue_manifold(residual_fn, u0, p0, alg)
    else
        throw(ArgumentError("Unsupported invariant set representation for continuation."))
    end
end

end # module Continuation