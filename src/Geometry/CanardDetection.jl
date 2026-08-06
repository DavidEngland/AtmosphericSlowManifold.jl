@enum SingularType FOLDED_NODE FOLDED_SADDLE FOLDED_FOCUS

struct FoldedSingularity
    point::ManifoldPoint
    classification::Symbol
    singular_type::Union{Nothing, SingularType}
    eigvals::Vector{ComplexF64}
end

struct CanardSegment <: AbstractInvariantSet
    singular_point::Vector{Float64}
    classification::SingularType
    canard_trajectory::Matrix{Float64}
    eigenvalues::Vector{ComplexF64}
end

function classify_folded_singularity(J::AbstractMatrix{<:Real}; tol::Float64 = 1e-10)
    eigs = eigvals(Matrix{Float64}(J))
    if length(eigs) != 2
        return :indeterminate
    end

    tr = real(sum(eigs))
    detJ = real(prod(eigs))
    disc = tr^2 - 4detJ

    if detJ < -tol
        return :folded_saddle
    elseif detJ > tol && disc >= 0
        return :folded_node
    elseif detJ > tol && disc < 0
        return :folded_focus
    else
        return :degenerate
    end
end

function classify_singular_type(J::AbstractMatrix{<:Real}; tol::Float64 = 1e-10)
    cls = classify_folded_singularity(J; tol = tol)
    if cls == :folded_node
        return FOLDED_NODE
    elseif cls == :folded_saddle
        return FOLDED_SADDLE
    elseif cls == :folded_focus
        return FOLDED_FOCUS
    end
    return nothing
end

function detect_folded_singularity(vfield, p::ManifoldPoint; eps::Float64 = 1e-6)
    y = p.y
    n = length(y)
    f0 = Float64.(vfield(p.x, y))
    J = zeros(n, n)

    for j in 1:n
        yp = copy(y)
        ym = copy(y)
        yp[j] += eps
        ym[j] -= eps
        fp = Float64.(vfield(p.x, yp))
        fm = Float64.(vfield(p.x, ym))
        J[:, j] = (fp .- fm) ./ (2eps)
    end

    cls = classify_folded_singularity(J)
    st = classify_singular_type(J)
    return FoldedSingularity(p, cls, st, ComplexF64.(eigvals(J)))
end

function build_canard_segment(
    folded::FoldedSingularity,
    trajectory::Matrix{Float64},
)
    folded.singular_type === nothing && throw(ArgumentError("Cannot build CanardSegment from degenerate or indeterminate folded singularity."))
    return CanardSegment(vcat(folded.point.x, folded.point.y), folded.singular_type, trajectory, folded.eigvals)
end
