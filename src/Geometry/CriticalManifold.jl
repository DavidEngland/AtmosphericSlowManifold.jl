abstract type AbstractInvariantSet end

struct ManifoldPoint
    x::Vector{Float64}
    y::Vector{Float64}
    residual_norm::Float64
end

struct CriticalManifoldSurface <: AbstractInvariantSet
    points::Vector{ManifoldPoint}
    coordinates::Matrix{Float64}
    fast_nullclines::Vector{Vector{Float64}}
    stability_mask::Vector{Bool}
    dimension::Int
end

function CriticalManifoldSurface(points::Vector{ManifoldPoint}, fast_nullclines::Vector{Vector{Float64}}, stability_mask::Vector{Bool})
    npts = length(points)
    npts == length(fast_nullclines) || throw(ArgumentError("fast_nullclines length must match points length."))
    npts == length(stability_mask) || throw(ArgumentError("stability_mask length must match points length."))
    npts == 0 && return CriticalManifoldSurface(points, zeros(0, 0), fast_nullclines, stability_mask, 0)

    dim = length(points[1].x) + length(points[1].y)
    coords = zeros(Float64, dim, npts)
    for (i, p) in enumerate(points)
        coords[:, i] = vcat(p.x, p.y)
    end
    return CriticalManifoldSurface(points, coords, fast_nullclines, stability_mask, dim)
end

function find_manifold_point(
    f_fast,
    x::AbstractVector{<:Real},
    y0::AbstractVector{<:Real};
    params = Dict{Symbol, Float64}(),
    tol::Float64 = 1e-8,
    maxiter::Int = 40,
)
    y = Float64.(y0)
    xv = Float64.(x)

    for _ in 1:maxiter
        r = Float64.(f_fast(xv, y, 0.0; params = params))
        rn = norm(r)
        if rn <= tol
            return ManifoldPoint(xv, y, rn)
        end

        J = finite_difference_jacobian_y(f_fast, xv, y; params = params)
        delta = J \ r
        y .-= delta

        if norm(delta) <= tol
            r2 = Float64.(f_fast(xv, y, 0.0; params = params))
            return ManifoldPoint(xv, y, norm(r2))
        end
    end

    r = Float64.(f_fast(xv, y, 0.0; params = params))
    return ManifoldPoint(xv, y, norm(r))
end

function solve_critical_surface(
    f_fast,
    x_path::Vector{<:AbstractVector{<:Real}},
    y_seed::AbstractVector{<:Real};
    params = Dict{Symbol, Float64}(),
    tol::Float64 = 1e-8,
    maxiter::Int = 40,
)
    points = ManifoldPoint[]
    fast_nullclines = Vector{Vector{Float64}}()
    stability_mask = Bool[]
    y = Float64.(y_seed)

    for x in x_path
        p = find_manifold_point(f_fast, x, y; params = params, tol = tol, maxiter = maxiter)
        push!(points, p)
        r = Float64.(f_fast(Float64.(p.x), Float64.(p.y), 0.0; params = params))
        push!(fast_nullclines, r)

        J = finite_difference_jacobian_y(f_fast, p.x, p.y; params = params)
        eigs = eigvals(Matrix{Float64}(J))
        push!(stability_mask, all(real.(eigs) .< 0.0))

        y = p.y
    end

    return CriticalManifoldSurface(points, fast_nullclines, stability_mask)
end
