module Jacobians

using LinearAlgebra
using Symbolics
using ModelingToolkit

export JacobianModel, JacobianCache
export evaluate_jacobian, evaluate_det, evaluate_adjugate
export compute_fast_jacobian, compute_adjugate, evaluate_tangent_space

struct JacobianModel
    fast_vars::Vector{Num}
    slow_vars::Vector{Num}
    params::Vector{Num}
    J_symbolic::Matrix{Num}
    det_symbolic::Num
    adj_symbolic::Matrix{Num}
    J_function::Function
    det_function::Function
    adj_function::Function
end

const JacobianCache = JacobianModel

function _compile_symbolic_kernel(expr, args...)
    compiled = Symbolics.build_function(expr, args...; expression = Val(false))
    kernel = compiled isa Tuple ? compiled[1] : compiled
    kernel isa Expr && (kernel = eval(kernel))
    return kernel
end

function JacobianModel(
    f_fast::Vector{Num},
    fast_vars::Vector{Num},
    slow_vars::Vector{Num};
    params::Vector{Num} = Num[],
)
    J_sym = Symbolics.jacobian(f_fast, fast_vars)
    det_sym = Symbolics.det(J_sym)
    adj_sym = compute_adjugate(J_sym)

    J_fn = _compile_symbolic_kernel(J_sym, fast_vars, slow_vars, params)
    det_fn = _compile_symbolic_kernel(det_sym, fast_vars, slow_vars, params)
    adj_fn = _compile_symbolic_kernel(adj_sym, fast_vars, slow_vars, params)

    return JacobianModel(fast_vars, slow_vars, params, J_sym, det_sym, adj_sym, J_fn, det_fn, adj_fn)
end

function compute_adjugate(J::AbstractMatrix{<:Number})
    n, m = size(J)
    n == m || throw(ArgumentError("Adjugate requires a square matrix."))

    if n == 2
        return [J[2, 2] -J[1, 2]; -J[2, 1] J[1, 1]]
    elseif n == 3
        cof = Matrix{eltype(J)}(undef, 3, 3)
        cof[1, 1] = J[2, 2] * J[3, 3] - J[2, 3] * J[3, 2]
        cof[1, 2] = -(J[2, 1] * J[3, 3] - J[2, 3] * J[3, 1])
        cof[1, 3] = J[2, 1] * J[3, 2] - J[2, 2] * J[3, 1]
        cof[2, 1] = -(J[1, 2] * J[3, 3] - J[1, 3] * J[3, 2])
        cof[2, 2] = J[1, 1] * J[3, 3] - J[1, 3] * J[3, 1]
        cof[2, 3] = -(J[1, 1] * J[3, 2] - J[1, 2] * J[3, 1])
        cof[3, 1] = J[1, 2] * J[2, 3] - J[1, 3] * J[2, 2]
        cof[3, 2] = -(J[1, 1] * J[2, 3] - J[1, 3] * J[2, 1])
        cof[3, 3] = J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]
        return transpose(cof)
    elseif eltype(J) <: Num
        throw(ArgumentError("Symbolic adjugate currently supports only 2x2 and 3x3 matrices."))
    else
        return det(J) * inv(J)
    end
end

function compute_fast_jacobian(
    f_fast,
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real};
    eps::Float64 = 1e-6,
    params = Dict{Symbol, Float64}(),
)
    y0 = Float64.(y)
    f0 = Float64.(f_fast(Float64.(x), y0, 0.0; params = params))
    m = length(f0)
    n = length(y0)
    J = zeros(m, n)

    for j in 1:n
        yp = copy(y0)
        ym = copy(y0)
        yp[j] += eps
        ym[j] -= eps
        fp = Float64.(f_fast(Float64.(x), yp, 0.0; params = params))
        fm = Float64.(f_fast(Float64.(x), ym, 0.0; params = params))
        J[:, j] = (fp .- fm) ./ (2eps)
    end

    return J
end

function _vectorize_values(vars::Vector{Num}, substitutions::Dict{Num, <:Real})
    vals = Float64[]
    for v in vars
        haskey(substitutions, v) || throw(ArgumentError("Missing substitution for symbolic variable $(v)."))
        push!(vals, float(substitutions[v]))
    end
    return vals
end

function evaluate_jacobian(model::JacobianModel, u_fast::AbstractVector, u_slow::AbstractVector, p::AbstractVector)
    return Matrix{Float64}(model.J_function(u_fast, u_slow, p))
end

function evaluate_det(model::JacobianModel, u_fast::AbstractVector, u_slow::AbstractVector, p::AbstractVector)
    return float(model.det_function(u_fast, u_slow, p))
end

function evaluate_adjugate(model::JacobianModel, u_fast::AbstractVector, u_slow::AbstractVector, p::AbstractVector)
    return Matrix{Float64}(model.adj_function(u_fast, u_slow, p))
end

function compute_fast_jacobian(model::JacobianModel, substitutions::Dict{Num, <:Real})
    u_fast = _vectorize_values(model.fast_vars, substitutions)
    u_slow = _vectorize_values(model.slow_vars, substitutions)
    p = _vectorize_values(model.params, substitutions)
    return evaluate_jacobian(model, u_fast, u_slow, p)
end

function evaluate_tangent_space(J::AbstractMatrix{<:Real}; atol::Float64 = 1e-10)
    F = svd(Matrix{Float64}(J))
    if isempty(F.S)
        return zeros(size(J, 2), 0)
    end
    maxs = maximum(F.S)
    thresh = max(atol, atol * maxs)
    idx = findall(s -> s <= thresh, F.S)
    if isempty(idx)
        return zeros(size(J, 2), 0)
    end
    return F.V[:, idx]
end

end # module Jacobians