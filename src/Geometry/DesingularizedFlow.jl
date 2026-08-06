function slow_flow_vector(f_fast, g_slow, x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; params = Dict{Symbol, Float64}())
    J = finite_difference_jacobian_y(f_fast, x, y; params = params)
    g = Float64.(g_slow(Float64.(x), Float64.(y), 0.0; params = params))
    return J \ g
end

function desingularized_vector_field(f_fast, g_slow, x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; params = Dict{Symbol, Float64}())
    J = finite_difference_jacobian_y(f_fast, x, y; params = params)
    g = Float64.(g_slow(Float64.(x), Float64.(y), 0.0; params = params))
    return compute_adjugate(J) * g
end

function find_desingularized_singular_points(
    surface::CriticalManifoldSurface,
    f_fast,
    g_slow;
    params = Dict{Symbol, Float64}(),
    flow_tol::Float64 = 1e-6,
    det_tol::Float64 = 1e-6,
)
    out = ManifoldPoint[]
    for p in surface.points
        d = abs(fold_indicator(f_fast, p.x, p.y; params = params))
        v = desingularized_vector_field(f_fast, g_slow, p.x, p.y; params = params)
        if d <= det_tol && norm(v) <= flow_tol
            push!(out, p)
        end
    end
    return out
end
