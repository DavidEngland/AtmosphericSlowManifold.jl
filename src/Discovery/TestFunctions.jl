abstract type AbstractTestFunctionFamily end

struct GegenbauerFamily <: AbstractTestFunctionFamily
    lambda::Float64
    max_mode::Int
end

struct BSplineFamily <: AbstractTestFunctionFamily
    order::Int
    num_knots::Int
end

function _tf_gegenbauerC(n::Int, lambda::Float64, x::Float64)
    n == 0 && return 1.0
    n == 1 && return 2.0 * lambda * x

    c_nm2 = 1.0
    c_nm1 = 2.0 * lambda * x
    for k in 2:n
        c_n = (2.0 * (k + lambda - 1.0) * x * c_nm1 - (k + 2.0 * lambda - 2.0) * c_nm2) / k
        c_nm2 = c_nm1
        c_nm1 = c_n
    end
    return c_nm1
end

function _ref_coordinate(x::Float64, x0::Float64, x1::Float64)
    x1 > x0 || throw(ArgumentError("Reference interval must be non-degenerate."))
    return 2.0 * (x - x0) / (x1 - x0) - 1.0
end

function evaluate_test_function(family::GegenbauerFamily, mode::Int, z::Float64, z0::Float64, H::Float64)
    mode >= 0 || throw(ArgumentError("mode must be non-negative"))
    z_ref = _ref_coordinate(z, z0, H)
    return _tf_gegenbauerC(mode, family.lambda, z_ref) * (1.0 - z_ref^2)^2
end

function evaluate_dt_test_function(::GegenbauerFamily, mode::Int, t::Float64, t0::Float64, T::Float64)
    T > t0 || throw(ArgumentError("Time interval must be non-degenerate."))
    mode >= 1 || return 0.0
    τ = (t - t0) / (T - t0)
    return (mode - 1) == 0 ? 0.0 : (mode - 1) * τ^(mode - 2) / (T - t0)
end

function evaluate_dz2_test_function(family::GegenbauerFamily, mode::Int, z::Float64, z0::Float64, H::Float64)
    h = max(1e-6, 1e-6 * (H - z0))
    z1 = clamp(z - h, z0, H)
    z2 = clamp(z + h, z0, H)
    f1 = evaluate_test_function(family, mode, z1, z0, H)
    f0 = evaluate_test_function(family, mode, z, z0, H)
    f2 = evaluate_test_function(family, mode, z2, z0, H)
    denom = max((z2 - z)^2, 1e-12)
    return (f2 - 2f0 + f1) / denom
end

function evaluate_test_function(family::BSplineFamily, mode::Int, z::Float64, z0::Float64, H::Float64)
    family.num_knots > 1 || throw(ArgumentError("num_knots must exceed 1"))
    mode >= 1 || throw(ArgumentError("BSpline mode is 1-based"))

    ξ = _ref_coordinate(z, z0, H)
    knots = range(-1.0, 1.0; length = family.num_knots)
    kidx = clamp(mode, 1, length(knots))
    center = knots[kidx]
    width = 2.0 / max(1, family.num_knots - 1)
    r = abs(ξ - center) / max(width, 1e-8)
    return r < 1.0 ? (1.0 - r)^family.order : 0.0
end

function evaluate_dt_test_function(::BSplineFamily, mode::Int, t::Float64, t0::Float64, T::Float64)
    T > t0 || throw(ArgumentError("Time interval must be non-degenerate."))
    mode >= 1 || return 0.0
    τ = (t - t0) / (T - t0)
    return (mode - 1) == 0 ? 0.0 : (mode - 1) * τ^(mode - 2) / (T - t0)
end

function evaluate_dz2_test_function(family::BSplineFamily, mode::Int, z::Float64, z0::Float64, H::Float64)
    h = max(1e-6, 1e-6 * (H - z0))
    z1 = clamp(z - h, z0, H)
    z2 = clamp(z + h, z0, H)
    f1 = evaluate_test_function(family, mode, z1, z0, H)
    f0 = evaluate_test_function(family, mode, z, z0, H)
    f2 = evaluate_test_function(family, mode, z2, z0, H)
    denom = max((z2 - z)^2, 1e-12)
    return (f2 - 2f0 + f1) / denom
end
