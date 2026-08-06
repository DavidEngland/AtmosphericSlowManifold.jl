struct FoldCurve <: AbstractInvariantSet
    points::Matrix{Float64}
    tangents::Matrix{Float64}
    transversality::Vector{Float64}
end

function fold_indicator(f_fast, x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; params = Dict{Symbol, Float64}())
    J = finite_difference_jacobian_y(f_fast, x, y; params = params)
    size(J, 1) == size(J, 2) || throw(ArgumentError("Fold indicator requires square D_y f."))
    return det(J)
end

function _newton_scalar_correct(residual_fn, u0::Vector{Float64}, p::Float64; tol::Float64 = 1e-7, maxiter::Int = 20)
    u = copy(u0)
    for _ in 1:maxiter
        r = residual_fn(u, p)
        norm(r) <= tol && return u, true

        # Local finite-difference Jacobian wrt state vector.
        n = length(u)
        m = length(r)
        m == n || return u, false
        J = zeros(m, n)
        h = 1e-6
        for j in 1:n
            up = copy(u)
            um = copy(u)
            up[j] += h
            um[j] -= h
            rp = residual_fn(up, p)
            rm = residual_fn(um, p)
            J[:, j] = (rp .- rm) ./ (2h)
        end

        du = J \ r
        u .-= du
        norm(du) <= tol && return u, true
    end
    return u, false
end

"""
Track a fold line for a parameterized residual function.

`residual_fn(u, p)` should return the fast residual at state `u` and parameter `p`.
The fold condition is appended as `det(D_u residual_fn) = 0`.
"""
function track_fold_line(
    residual_fn,
    u_init::Vector{Float64},
    p_range::Tuple{Float64, Float64};
    n_steps::Int = 32,
    tol::Float64 = 1e-6,
)
    pmin, pmax = p_range
    pvals = collect(range(pmin, pmax; length = max(2, n_steps)))

    point_list = Vector{Vector{Float64}}()
    tangent_list = Vector{Vector{Float64}}()
    trans = Float64[]
    u = copy(u_init)

    for p in pvals
        augmented(uvec, pval) = begin
            r = Float64.(residual_fn(uvec, pval))
            n = length(uvec)
            J = zeros(n, n)
            h = 1e-6
            for j in 1:n
                up = copy(uvec)
                um = copy(uvec)
                up[j] += h
                um[j] -= h
                rp = Float64.(residual_fn(up, pval))
                rm = Float64.(residual_fn(um, pval))
                J[:, j] = (rp .- rm) ./ (2h)
            end
            return [r; det(J)]
        end

        u_aug = [u; 0.0]
        corrected, ok = _newton_scalar_correct(augmented, u_aug, p; tol = tol)
        if ok
            u = corrected[1:length(u)]
            push!(point_list, [u; p])

            # Tangent proxy from local Jacobian nullspace direction.
            n = length(u)
            J = zeros(n, n)
            h_local = 1e-6
            for j in 1:n
                up = copy(u)
                um = copy(u)
                up[j] += h_local
                um[j] -= h_local
                rp = Float64.(residual_fn(up, p))
                rm = Float64.(residual_fn(um, p))
                J[:, j] = (rp .- rm) ./ (2h_local)
            end
            S = svd(J)
            t = S.V[:, end]
            push!(tangent_list, Float64.(t))

            # Numerical transversality proxy: d/dp det(D_u f).
            h = 1e-5
            det_at(pp) = begin
                n = length(u)
                J = zeros(n, n)
                for j in 1:n
                    up = copy(u)
                    um = copy(u)
                    up[j] += 1e-6
                    um[j] -= 1e-6
                    rp = Float64.(residual_fn(up, pp))
                    rm = Float64.(residual_fn(um, pp))
                    J[:, j] = (rp .- rm) ./ (2e-6)
                end
                return det(J)
            end
            τ = (det_at(p + h) - det_at(p - h)) / (2h)
            push!(trans, τ)
        end
    end

    if isempty(point_list)
        return FoldCurve(zeros(length(u_init) + 1, 0), zeros(length(u_init), 0), Float64[])
    end

    pts = reduce(hcat, point_list)
    tans = reduce(hcat, tangent_list)
    return FoldCurve(pts, tans, trans)
end

function track_fold_curve(
    f_fast,
    x_path::Vector{<:AbstractVector{<:Real}},
    y_seed::AbstractVector{<:Real};
    params = Dict{Symbol, Float64}(),
    tol::Float64 = 1e-6,
)
    surface = solve_critical_surface(f_fast, x_path, y_seed; params = params)
    fold_points = ManifoldPoint[]

    for p in surface.points
        σ = fold_indicator(f_fast, p.x, p.y; params = params)
        if abs(σ) <= tol
            push!(fold_points, p)
        end
    end

    return fold_points
end
