using LinearAlgebra
using Statistics

struct IdentifiabilityReport{T<:Real}
    gram_matrix::Matrix{T}
    mutual_coherence::T
    condition_number::T
    singular_values::Vector{T}
    vif::Vector{T}
    parameter_covariance::Matrix{T}
end

function compute_gram_matrix(G::AbstractMatrix{<:Real})
    return Matrix{Float64}(G' * G)
end

function compute_mutual_coherence(G::AbstractMatrix{<:Real}; eps_floor::Float64 = 1e-12)
    X = Matrix{Float64}(G)
    n = size(X, 2)
    n <= 1 && return 0.0

    norms = vec(sqrt.(sum(abs2, X; dims = 1)))
    Cmax = 0.0
    for i in 1:n
        for j in (i + 1):n
            denom = max(norms[i] * norms[j], eps_floor)
            cij = abs(dot(X[:, i], X[:, j])) / denom
            if cij > Cmax
                Cmax = cij
            end
        end
    end
    return Cmax
end

function compute_condition_number(G::AbstractMatrix{<:Real}; eps_floor::Float64 = 1e-12)
    s = svdvals(Matrix{Float64}(G))
    isempty(s) && return Inf
    return maximum(s) / max(minimum(s), eps_floor)
end

function compute_vif(G::AbstractMatrix{<:Real}; eps_floor::Float64 = 1e-12)
    X = Matrix{Float64}(G)
    p = size(X, 2)
    p == 0 && return Float64[]
    p == 1 && return [1.0]

    out = zeros(Float64, p)
    for j in 1:p
        y = X[:, j]
        cols = [k for k in 1:p if k != j]
        Xo = X[:, cols]

        beta = Xo \ y
        yhat = Xo * beta
        rss = sum(abs2, y .- yhat)
        tss = sum(abs2, y .- mean(y))

        if tss <= eps_floor
            out[j] = Inf
        else
            r2 = clamp(1.0 - rss / tss, 0.0, 1.0 - eps_floor)
            out[j] = 1.0 / max(1.0 - r2, eps_floor)
        end
    end
    return out
end

function compute_parameter_covariance(
    G::AbstractMatrix{<:Real};
    sigma2::Float64 = 1.0,
    ridge::Float64 = 0.0,
)
    sigma2 >= 0 || throw(ArgumentError("sigma2 must be non-negative."))
    ridge >= 0 || throw(ArgumentError("ridge must be non-negative."))

    X = Matrix{Float64}(G)
    p = size(X, 2)
    H = X' * X + ridge * I(p)
    return sigma2 .* pinv(H)
end

function analyze_identifiability(
    G::AbstractMatrix{<:Real};
    sigma2::Float64 = 1.0,
    ridge::Float64 = 0.0,
)
    gram = compute_gram_matrix(G)
    svals = svdvals(Matrix{Float64}(G))
    return IdentifiabilityReport(
        gram,
        compute_mutual_coherence(G),
        compute_condition_number(G),
        Float64.(svals),
        compute_vif(G),
        compute_parameter_covariance(G; sigma2 = sigma2, ridge = ridge),
    )
end

function prune_by_mutual_coherence(
    G::AbstractMatrix{<:Real};
    threshold::Float64 = 0.85,
)
    0.0 <= threshold <= 1.0 || throw(ArgumentError("threshold must be in [0,1]."))

    X = Matrix{Float64}(G)
    p = size(X, 2)
    if p <= 1
        return (kept = collect(1:p), dropped = Int[])
    end

    norms = vec(sqrt.(sum(abs2, X; dims = 1)))
    kept = Int[]
    dropped = Int[]

    for j in 1:p
        if norms[j] <= eps(Float64)
            push!(dropped, j)
            continue
        end

        too_correlated = false
        for k in kept
            denom = max(norms[j] * norms[k], eps(Float64))
            cij = abs(dot(X[:, j], X[:, k])) / denom
            if cij > threshold
                too_correlated = true
                break
            end
        end

        if too_correlated
            push!(dropped, j)
        else
            push!(kept, j)
        end
    end

    return (kept = kept, dropped = dropped)
end
