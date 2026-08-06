struct HyperbolicityReport
    eigvals::Vector{ComplexF64}
    min_abs_real::Float64
    spectral_gap::Float64
    epsilon0::Float64
    is_hyperbolic::Bool
end

function fenichel_metrics(J::AbstractMatrix{<:Real}; gap_tol::Float64 = 1e-6)
    eigs = ComplexF64.(eigvals(Matrix{Float64}(J)))
    reals = abs.(real.(eigs))
    sorted = sort(reals)
    min_abs_real = sorted[1]
    gap = length(sorted) >= 2 ? sorted[2] - sorted[1] : sorted[1]
    eps0 = max(0.0, min_abs_real)
    is_hyp = min_abs_real > gap_tol
    return HyperbolicityReport(eigs, min_abs_real, gap, eps0, is_hyp)
end

function fenichel_metrics(f_fast, x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; params = Dict{Symbol, Float64}(), gap_tol::Float64 = 1e-6)
    J = finite_difference_jacobian_y(f_fast, x, y; params = params)
    return fenichel_metrics(J; gap_tol = gap_tol)
end

function fenichel_metrics(
    model::JacobianModel,
    u_fast::AbstractVector{<:Real},
    u_slow::AbstractVector{<:Real},
    p::AbstractVector{<:Real};
    gap_tol::Float64 = 1e-6,
)
    J = evaluate_jacobian(model, Float64.(u_fast), Float64.(u_slow), Float64.(p))
    return fenichel_metrics(J; gap_tol = gap_tol)
end

function hyperbolicity_profile(
    surface::CriticalManifoldSurface,
    f_fast;
    params = Dict{Symbol, Float64}(),
    gap_tol::Float64 = 1e-6,
)
    return [fenichel_metrics(f_fast, p.x, p.y; params = params, gap_tol = gap_tol) for p in surface.points]
end
