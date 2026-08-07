"""
    smooth_max(x; eps=1e-3)

Smooth approximation of `max(x, 0)` with C-infinity regularity.
"""
@inline function smooth_max(x::TX; eps::TE = 1e-3) where {TX<:Real, TE<:Real}
    x1, e1 = promote(x, eps)
    return (x1 + sqrt(x1 * x1 + e1 * e1)) / 2
end

"""
    smooth_min(x; eps=1e-3)

Smooth approximation of `min(x, 0)` with C-infinity regularity.
"""
@inline function smooth_min(x::TX; eps::TE = 1e-3) where {TX<:Real, TE<:Real}
    x1, e1 = promote(x, eps)
    return (x1 - sqrt(x1 * x1 + e1 * e1)) / 2
end

"""
    smooth_floor(x, floor_val; eps=1e-3)

Smooth lower-bound clamp approximation of `max(x, floor_val)`.
"""
@inline function smooth_floor(x::TX, floor_val::TF; eps::TE = 1e-3) where {TX<:Real, TF<:Real, TE<:Real}
    x1, f1, e1 = promote(x, floor_val, eps)
    return f1 + smooth_max(x1 - f1; eps = e1)
end