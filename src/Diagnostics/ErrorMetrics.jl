module ErrorMetrics

using LinearAlgebra
using Statistics

export rmse, mae, bias, r2, nrmse, skill_score, correlation,
       normalized_bias, relative_l2_error, closure_residual

@inline function _assert_same_length(y_pred, y_true)
    length(y_pred) == length(y_true) || throw(ArgumentError("Array dimensions must match"))
end

"""
    rmse(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where {T<:Real}

Computes Root Mean Square Error.
"""
function rmse(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where {T<:Real}
    _assert_same_length(y_pred, y_true)
    n = length(y_true)
    n == 0 && return zero(Float64)

    acc = 0.0
    @inbounds for i in eachindex(y_true, y_pred)
        d = Float64(y_pred[i]) - Float64(y_true[i])
        acc += d * d
    end
    return sqrt(acc / n)
end

"""
    mae(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where {T<:Real}

Computes Mean Absolute Error.
"""
function mae(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where {T<:Real}
    _assert_same_length(y_pred, y_true)
    n = length(y_true)
    n == 0 && return zero(Float64)

    acc = 0.0
    @inbounds for i in eachindex(y_true, y_pred)
        acc += abs(Float64(y_pred[i]) - Float64(y_true[i]))
    end
    return acc / n
end

"""
    bias(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where {T<:Real}

Computes mean systematic bias.
"""
function bias(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where {T<:Real}
    _assert_same_length(y_pred, y_true)
    n = length(y_true)
    n == 0 && return zero(Float64)

    acc = 0.0
    @inbounds for i in eachindex(y_true, y_pred)
        acc += Float64(y_pred[i]) - Float64(y_true[i])
    end
    return acc / n
end

"""
    r2(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where {T<:Real}

Computes coefficient of determination.
"""
function r2(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where {T<:Real}
    _assert_same_length(y_pred, y_true)
    n = length(y_true)
    n == 0 && return one(Float64)

    ybar = mean(y_true)
    ss_res = 0.0
    ss_tot = 0.0
    @inbounds for i in eachindex(y_true, y_pred)
        yt = Float64(y_true[i])
        yp = Float64(y_pred[i])
        d1 = yt - yp
        d2 = yt - Float64(ybar)
        ss_res += d1 * d1
        ss_tot += d2 * d2
    end
    return iszero(ss_tot) ? 1.0 : 1.0 - (ss_res / ss_tot)
end

"""
    nrmse(y_pred::AbstractArray{T}, y_true::AbstractArray{T}; norm::Symbol=:std) where {T<:Real}

Computes normalized RMSE using `:std`, `:range`, or `:mean` scale.
"""
function nrmse(y_pred::AbstractArray{T}, y_true::AbstractArray{T}; norm::Symbol = :std) where {T<:Real}
    val_rmse = rmse(y_pred, y_true)
    scale = if norm == :std
        std(y_true)
    elseif norm == :range
        maximum(y_true) - minimum(y_true)
    elseif norm == :mean
        abs(mean(y_true))
    else
        throw(ArgumentError("Unsupported normalization option: $(norm)"))
    end
    scale_f = Float64(scale)
    return iszero(scale_f) ? val_rmse : val_rmse / scale_f
end

"""
    skill_score(y_pred::AbstractArray{T}, y_true::AbstractArray{T}, y_ref::AbstractArray{T}) where {T<:Real}

Computes Murphy skill score relative to a reference trajectory.
"""
function skill_score(y_pred::AbstractArray{T}, y_true::AbstractArray{T}, y_ref::AbstractArray{T}) where {T<:Real}
    _assert_same_length(y_pred, y_true)
    _assert_same_length(y_ref, y_true)
    n = length(y_true)
    n == 0 && return zero(Float64)

    mse_pred = 0.0
    mse_ref = 0.0
    @inbounds for i in eachindex(y_true, y_pred, y_ref)
        yt = Float64(y_true[i])
        dp = Float64(y_pred[i]) - yt
        dr = Float64(y_ref[i]) - yt
        mse_pred += dp * dp
        mse_ref += dr * dr
    end
    mse_pred /= n
    mse_ref /= n
    return iszero(mse_ref) ? 0.0 : 1.0 - (mse_pred / mse_ref)
end

"""
    correlation(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where {T<:Real}

Computes Pearson correlation coefficient between prediction and truth.
"""
function correlation(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where {T<:Real}
    _assert_same_length(y_pred, y_true)
    n = length(y_true)
    n == 0 && return zero(Float64)

    μp = mean(y_pred)
    μt = mean(y_true)
    num = 0.0
    den_p = 0.0
    den_t = 0.0

    @inbounds for i in eachindex(y_true, y_pred)
        dp = Float64(y_pred[i]) - Float64(μp)
        dt = Float64(y_true[i]) - Float64(μt)
        num += dp * dt
        den_p += dp * dp
        den_t += dt * dt
    end
    den = sqrt(den_p * den_t)
    return iszero(den) ? 0.0 : num / den
end

"""
    normalized_bias(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where {T<:Real}

Computes normalized bias: `bias / abs(mean(y_true))`.
"""
function normalized_bias(y_pred::AbstractArray{T}, y_true::AbstractArray{T}) where {T<:Real}
    b = bias(y_pred, y_true)
    μ = abs(Float64(mean(y_true)))
    return iszero(μ) ? b : b / μ
end

"""
    relative_l2_error(y_pred::AbstractVector{T}, y_true::AbstractVector{T}) where {T<:Real}

Computes relative L2 error norm.
"""
function relative_l2_error(y_pred::AbstractVector{T}, y_true::AbstractVector{T}) where {T<:Real}
    _assert_same_length(y_pred, y_true)
    num = 0.0
    den = 0.0
    @inbounds for i in eachindex(y_true, y_pred)
        dy = Float64(y_pred[i]) - Float64(y_true[i])
        yt = Float64(y_true[i])
        num += dy * dy
        den += yt * yt
    end
    return iszero(den) ? sqrt(num) : sqrt(num / den)
end

"""
    closure_residual(phi_obs::AbstractVector{T}, phi_model::AbstractVector{T}) where {T<:Real}

Computes pointwise absolute closure residuals.
"""
function closure_residual(phi_obs::AbstractVector{T}, phi_model::AbstractVector{T}) where {T<:Real}
    _assert_same_length(phi_obs, phi_model)
    out = similar(phi_obs, Float64)
    @inbounds for i in eachindex(phi_obs, phi_model)
        out[i] = abs(Float64(phi_obs[i]) - Float64(phi_model[i]))
    end
    return out
end

end # module ErrorMetrics
