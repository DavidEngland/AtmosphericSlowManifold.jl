module ManifoldMetrics

using LinearAlgebra

export transversality, fold_distance, slow_manifold_error, normal_hyperbolicity

"""
    transversality(v_fast, n_slow)

Computes normalized alignment between fast tendency and slow-manifold normal.
"""
function transversality(v_fast::AbstractVector{T}, n_slow::AbstractVector{T}) where {T<:Real}
    length(v_fast) == length(n_slow) || throw(ArgumentError("Dimensions must match"))
    num = 0.0
    nv = 0.0
    nn = 0.0
    @inbounds for i in eachindex(v_fast, n_slow)
        v = Float64(v_fast[i])
        n = Float64(n_slow[i])
        num += v * n
        nv += v * v
        nn += n * n
    end
    den = sqrt(nv * nn)
    return iszero(den) ? 0.0 : abs(num) / den
end

"""
    fold_distance(state, fold_locus)

Computes distance from state vector to nearest fold-locus column.
"""
function fold_distance(state::AbstractVector{T}, fold_locus::AbstractMatrix{T}) where {T<:Real}
    length(state) == size(fold_locus, 1) || throw(ArgumentError("State and fold locus dimensions are incompatible."))
    m = size(fold_locus, 2)
    m == 0 && return typemax(Float64)

    min_dist2 = typemax(Float64)
    @inbounds for j in 1:m
        d2 = 0.0
        @views for i in eachindex(state)
            diff = Float64(state[i]) - Float64(fold_locus[i, j])
            d2 += diff * diff
        end
        if d2 < min_dist2
            min_dist2 = d2
        end
    end
    return sqrt(min_dist2)
end

"""
    slow_manifold_error(state, manifold_map)

Computes residual distance between state and mapped slow-manifold state.
"""
function slow_manifold_error(state::AbstractVector{T}, manifold_map::Function) where {T<:Real}
    mapped = manifold_map(state)
    length(mapped) == length(state) || throw(ArgumentError("Mapped manifold state dimension mismatch."))
    acc = 0.0
    @inbounds for i in eachindex(state)
        d = Float64(state[i]) - Float64(mapped[i])
        acc += d * d
    end
    return sqrt(acc)
end

"""
    normal_hyperbolicity(jacobian, n_fast)

Computes spectral-gap ratio between fast and slow eigenvalue magnitudes.
"""
function normal_hyperbolicity(jacobian::AbstractMatrix{T}, n_fast::Int) where {T<:Real}
    vals = eigen(jacobian).values
    rates = sort(abs.(real.(vals)); rev = true)
    (n_fast >= 1 && n_fast < length(rates)) || throw(ArgumentError("n_fast must satisfy 1 <= n_fast < n_eigs"))
    fast_rate = rates[n_fast]
    slow_rate = rates[n_fast + 1]
    return iszero(slow_rate) ? typemax(Float64) : fast_rate / slow_rate
end

end # module ManifoldMetrics
