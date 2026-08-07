using ForwardDiff
using LinearAlgebra

struct FoldedSingularity{T<:Real}
    fast_state::Vector{T}
    slow_state::Vector{T}
    classification::Symbol
    eigenvalues::Vector{Complex{T}}
    det_fast_jacobian::T
    spectral_ratio::T
end

function _adjugate_matrix(A::AbstractMatrix{T}) where {T<:Number}
    n, m = size(A)
    n == m || throw(ArgumentError("Adjugate requires a square matrix."))

    if n == 1
        return reshape([one(T)], 1, 1)
    end

    adj = Matrix{T}(undef, n, n)
    for i in 1:n
        for j in 1:n
            rows = [r for r in 1:n if r != j]
            cols = [c for c in 1:n if c != i]
            minor = A[rows, cols]
            sign = isodd(i + j) ? -one(T) : one(T)
            adj[i, j] = sign * det(minor)
        end
    end
    return adj
end

function _classify_from_eigenvalues(eigs::AbstractVector{<:Complex}; tol::Float64 = 1e-8)
    nonzero = [λ for λ in eigs if abs(λ) > tol]
    if length(nonzero) < 2
        return :nonsingular, 0.0
    end

    sort!(nonzero; by = x -> -abs(x))
    λ1, λ2 = nonzero[1], nonzero[2]
    ratio = abs(λ1) / max(abs(λ2), tol)

    if abs(imag(λ1)) > tol || abs(imag(λ2)) > tol
        return :folded_focus, ratio
    end

    prod_re = real(λ1) * real(λ2)
    if prod_re < -tol
        return :folded_saddle, ratio
    elseif prod_re > tol
        return :folded_node, ratio
    end

    return :degenerate, ratio
end

function classify_folded_singularity(J::AbstractMatrix{<:Real}; tol::Float64 = 1e-8)
    eigs = ComplexF64.(eigvals(Matrix{Float64}(J)))
    cls, _ = _classify_from_eigenvalues(eigs; tol = tol)
    return cls
end

function _desingularized_rhs(
    y::AbstractVector,
    z::AbstractVector,
    g_fn::Function,
    f_fn::Function,
)
    J_y = ForwardDiff.jacobian(yv -> g_fn(yv, z), y)
    J_z = ForwardDiff.jacobian(zv -> g_fn(y, zv), z)
    size(J_y, 1) == size(J_y, 2) || throw(ArgumentError("D_y g must be square for adjugate formulation."))

    f_val = f_fn(y, z)
    length(f_val) == size(J_z, 2) || throw(ArgumentError("f(y,z) dimension must match number of slow variables in D_z g."))

    adj_J_y = _adjugate_matrix(Matrix(J_y))
    return -(adj_J_y * Matrix(J_z) * collect(f_val))
end

function classify_folded_singularity(
    y_star::AbstractVector{<:Real},
    z_star::AbstractVector{<:Real},
    g_fn::Function,
    f_fn::Function;
    tol::Float64 = 1e-8,
)
    y0 = Vector{Float64}(y_star)
    z0 = Vector{Float64}(z_star)
    ny = length(y0)

    x0 = vcat(y0, z0)
    full_rhs = function (x)
        y = x[1:ny]
        z = x[(ny + 1):end]
        ydot = _desingularized_rhs(y, z, g_fn, f_fn)
        zdot = collect(f_fn(y, z))
        return vcat(ydot, zdot)
    end

    J_desing = ForwardDiff.jacobian(full_rhs, x0)
    eigs = ComplexF64.(eigvals(Matrix{Float64}(J_desing)))
    cls, ratio = _classify_from_eigenvalues(eigs; tol = tol)

    J_y = ForwardDiff.jacobian(yv -> g_fn(yv, z0), y0)
    det_jy = det(Matrix{Float64}(J_y))

    return FoldedSingularity(y0, z0, cls, eigs, det_jy, ratio)
end

function detect_canard_trajectories(
    t::AbstractVector{<:Real},
    states::AbstractVector{<:AbstractVector{<:Real}},
    g_fn::Function;
    ny::Int,
    fold_threshold::Float64 = 1e-6,
    persistence_threshold::Float64 = 0.1,
)
    length(t) == length(states) || throw(ArgumentError("t and states must have the same length."))
    ny > 0 || throw(ArgumentError("ny must be positive."))

    segments = NamedTuple{(:t_start, :t_end, :duration, :start_index, :end_index)}[]
    in_repelling = false
    start_idx = 0

    for i in eachindex(t)
        u = states[i]
        length(u) > ny || throw(ArgumentError("Each state must contain both fast and slow variables."))
        y = Vector{Float64}(u[1:ny])
        z = Vector{Float64}(u[(ny + 1):end])
        J_y = ForwardDiff.jacobian(yv -> g_fn(yv, z), y)
        det_jy = det(Matrix{Float64}(J_y))

        on_repelling_sheet = det_jy < -abs(fold_threshold)
        if on_repelling_sheet && !in_repelling
            in_repelling = true
            start_idx = i
        elseif !on_repelling_sheet && in_repelling
            in_repelling = false
            duration = float(t[i] - t[start_idx])
            if duration >= persistence_threshold
                push!(segments, (
                    t_start = float(t[start_idx]),
                    t_end = float(t[i]),
                    duration = duration,
                    start_index = start_idx,
                    end_index = i,
                ))
            end
        end
    end

    if in_repelling
        last_idx = length(t)
        duration = float(t[last_idx] - t[start_idx])
        if duration >= persistence_threshold
            push!(segments, (
                t_start = float(t[start_idx]),
                t_end = float(t[last_idx]),
                duration = duration,
                start_index = start_idx,
                end_index = last_idx,
            ))
        end
    end

    return segments
end

function detect_canard_trajectories(sol, g_fn::Function; ny::Int, kwargs...)
    hasproperty(sol, :t) || throw(ArgumentError("Solution object must provide a t field."))
    hasproperty(sol, :u) || throw(ArgumentError("Solution object must provide a u field."))
    return detect_canard_trajectories(sol.t, sol.u, g_fn; ny = ny, kwargs...)
end
