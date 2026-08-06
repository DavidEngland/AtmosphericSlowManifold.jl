struct SpectralBLGalerkin <: AbstractDiscretization
    n_modes::Int
    lambda::Float64
    H::Float64
end

function SpectralBLGalerkin(; n_modes::Int = 12, lambda::Float64 = 0.75, H::Float64 = 3000.0)
    return SpectralBLGalerkin(n_modes, lambda, H)
end

function _gegenbauerC(n::Int, lambda::Float64, x::Float64)
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

function _x_to_z(x::Float64, H::Float64)
    return 0.5 * H * (x + 1.0)
end

function _quad_nodes_weights(n::Int)
    xs = range(-1.0, 1.0; length = n)
    dx = step(xs)
    ws = fill(dx, n)
    ws[1] *= 0.5
    ws[end] *= 0.5
    return collect(xs), ws
end

function _mass_matrix(n_modes::Int, lambda::Float64; n_quad::Int = 256)
    x, w = _quad_nodes_weights(n_quad)
    M = zeros(n_modes, n_modes)
    for i in 1:n_modes
        for j in 1:n_modes
            acc = 0.0
            ni = i - 1
            nj = j - 1
            for q in eachindex(x)
                xq = x[q]
                weight = (1.0 - xq^2)^(lambda - 0.5)
                acc += _gegenbauerC(ni, lambda, xq) * _gegenbauerC(nj, lambda, xq) * weight * w[q]
            end
            M[i, j] = acc
        end
    end
    return M
end

function _stiffness_matrix(n_modes::Int, lambda::Float64, H::Float64; n_quad::Int = 256)
    x, w = _quad_nodes_weights(n_quad)
    K = zeros(n_modes, n_modes)
    dzdx = 0.5 * H
    inv_dzdx = 1.0 / dzdx
    eps = 1e-6
    for i in 1:n_modes
        for j in 1:n_modes
            acc = 0.0
            ni = i - 1
            nj = j - 1
            for q in eachindex(x)
                xq = x[q]
                weight = (1.0 - xq^2)^(lambda - 0.5)

                # Numerical derivative in reference space then chain-rule to z.
                dCi_dx = (_gegenbauerC(ni, lambda, min(1.0, xq + eps)) - _gegenbauerC(ni, lambda, max(-1.0, xq - eps))) / (2 * eps)
                dCj_dx = (_gegenbauerC(nj, lambda, min(1.0, xq + eps)) - _gegenbauerC(nj, lambda, max(-1.0, xq - eps))) / (2 * eps)
                dCi_dz = dCi_dx * inv_dzdx
                dCj_dz = dCj_dx * inv_dzdx

                # Integrate in physical z-space: dz = dzdx * dx
                acc += dCi_dz * dCj_dz * weight * dzdx * w[q]
            end
            K[i, j] = acc
        end
    end
    return K
end

function dispatch_solve(
    disc::SpectralBLGalerkin,
    pde_sys::PDESystem,
    closure::AbstractClosure,
    tspan::Tuple{Float64, Float64};
    solver = Rodas5P(),
    kwargs...
)
    ModelingToolkit.@variables t
    ModelingToolkit.@variables a(t)[1:disc.n_modes]

    # Galerkin inner-product projection of a linear diffusive prototype operator:
    # M * da/dt = -K * a  =>  da/dt = L * a, where L = -(M \ K)
    M = _mass_matrix(disc.n_modes, disc.lambda)
    K = _stiffness_matrix(disc.n_modes, disc.lambda, disc.H)
    L = -(M \ K)

    eqs = Vector{Equation}(undef, disc.n_modes)
    for n in 1:disc.n_modes
        rhs = zero(a[1])
        for k in 1:disc.n_modes
            rhs += L[n, k] * a[k]
        end
        eqs[n] = ModelingToolkit.Differential(t)(a[n]) ~ rhs
    end

    ModelingToolkit.@named modal_ode = ModelingToolkit.ODESystem(eqs, t)
    sys = structural_simplify(modal_ode)

    prob = ODEProblem(sys, ones(disc.n_modes), tspan)
    return solve(prob, solver; kwargs...)
end
