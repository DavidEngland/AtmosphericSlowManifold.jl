struct SpectralBLGalerkin <: AbstractDiscretization
    n_modes::Int
    lambda::Float64
    H::Float64
    enable_nonlinear::Bool
    nonlinear_scale::Float64
    advection_response_scale::Float64
    diffusivity_response_scale::Float64
end

struct SpectralNonlinearTensors
    triple::Array{Float64, 3}
    advection::Array{Float64, 3}
    diffusion_flux::Array{Float64, 3}
end

"""
Decomposed modal RHS tendencies for spectral transport budget analysis.
"""
struct ModalBudgetDiagnostic
    linear::Vector{Float64}
    advection::Vector{Float64}
    diffusion::Vector{Float64}
    total::Vector{Float64}
end

function SpectralBLGalerkin(
    ;
    n_modes::Int = 12,
    lambda::Float64 = 0.75,
    H::Float64 = 3000.0,
    enable_nonlinear::Bool = true,
    nonlinear_scale::Float64 = 1.0,
    advection_response_scale::Float64 = 1.0,
    diffusivity_response_scale::Float64 = 1.0,
)
    return SpectralBLGalerkin(
        n_modes,
        lambda,
        H,
        enable_nonlinear,
        nonlinear_scale,
        advection_response_scale,
        diffusivity_response_scale,
    )
end

# Backward-compatible 3-arg constructor.
SpectralBLGalerkin(n_modes::Int, lambda::Float64, H::Float64) = SpectralBLGalerkin(n_modes, lambda, H, true, 1.0, 1.0, 1.0)

# Compatibility constructor with unused z0 argument used in notes/spec prototypes.
SpectralBLGalerkin(
    n_modes::Int,
    lambda::Float64,
    z0::Float64,
    H::Float64,
    enable_nonlinear::Bool,
    nonlinear_scale::Float64,
    advection_response_scale::Float64,
    diffusivity_response_scale::Float64,
) = SpectralBLGalerkin(
    n_modes,
    lambda,
    H,
    enable_nonlinear,
    nonlinear_scale,
    advection_response_scale,
    diffusivity_response_scale,
)

function _closure_scalar(expr::Num, state::ManifoldState)
    subs = Dict(
        state.eta1 => 0.0,
        state.eta2 => 0.0,
        state.eta3 => 0.0,
        state.r => 0.0,
        state.omega => 0.0,
        state.chi => 0.0,
        state.pi_g => 0.0,
        state.lambdamin => 1.0,
        state.u => 0.0,
        state.v => 0.0,
        state.theta => 0.0,
        state.q => 0.0,
        state.u_star => 0.3,
        state.z0 => 0.1,
    )
    v = Symbolics.value(Symbolics.substitute(expr, subs))
    if v isa Number
        return float(v)
    end
    parsed = tryparse(Float64, string(v))
    parsed === nothing && return 0.0
    return parsed
end

function _closure_nonlinear_strength(closure::AbstractClosure)
    state = ManifoldState()
    km = eddy_momentum(closure, state)
    kh = eddy_heat(closure, state)
    km0 = _closure_scalar(km, state)
    kh0 = _closure_scalar(kh, state)
    adv = abs(km0)
    diff = abs(kh0)
    return adv, diff
end

function _spectral_linear_rhs(L::AbstractMatrix{<:Real}, a::AbstractVector)
    n = length(a)
    rhs = zeros(Float64, n)
    for i in 1:n
        for j in 1:n
            rhs[i] += L[i, j] * a[j]
        end
    end
    return rhs
end

function _spectral_nonlinear_advection_rhs(
    tensors::SpectralNonlinearTensors,
    a::AbstractVector{<:Real},
    scale::Float64,
)
    n = length(a)
    rhs = zeros(Float64, n)
    for k in 1:n, i in 1:n, j in 1:n
        rhs[k] += scale * tensors.advection[k, i, j] * a[i] * a[j]
    end
    return rhs
end

function _spectral_nonlinear_diffusion_rhs(
    tensors::SpectralNonlinearTensors,
    a::AbstractVector{<:Real},
    scale::Float64,
)
    n = length(a)
    rhs = zeros(Float64, n)
    for k in 1:n, i in 1:n, j in 1:n
        rhs[k] += scale * tensors.diffusion_flux[k, i, j] * a[i] * a[j]
    end
    return rhs
end

function _spectral_nonlinear_diffusion_rhs(
    tensors::SpectralNonlinearTensors,
    u_hat::AbstractVector{<:Real},
    K_hat::AbstractVector{<:Real},
)
    n = length(u_hat)
    m = min(length(K_hat), n)
    k_eff = m == 0 ? 0.0 : sum(abs.(Float64.(K_hat[1:m]))) / m
    return _spectral_nonlinear_diffusion_rhs(tensors, Float64.(u_hat), k_eff)
end

function _spectral_nonlinear_rhs(
    tensors::SpectralNonlinearTensors,
    a::AbstractVector{<:Real},
    adv_scale::Float64,
    diff_scale::Float64,
)
    return _spectral_nonlinear_advection_rhs(tensors, a, adv_scale) .+
           _spectral_nonlinear_diffusion_rhs(tensors, a, diff_scale)
end

function _linear_modal_operator(disc::SpectralBLGalerkin)
    M = _mass_matrix(disc.n_modes, disc.lambda)
    K = _stiffness_matrix(disc.n_modes, disc.lambda, disc.H)
    return -(M \ K)
end

"""
    evaluate_modal_budget(u_hat, K_hat, disc, tensors)

Computes explicit modal budget terms for diagnostic logging and explainability.
"""
function evaluate_modal_budget(
    u_hat::Vector{Float64},
    K_hat::Vector{Float64},
    disc::SpectralBLGalerkin,
    tensors::SpectralNonlinearTensors,
)
    length(u_hat) == disc.n_modes || throw(ArgumentError("u_hat length must match disc.n_modes."))

    Lop = _linear_modal_operator(disc)
    f_lin = _spectral_linear_rhs(Lop, u_hat)
    f_adv = _spectral_nonlinear_advection_rhs(tensors, u_hat, 1.0)
    f_diff = _spectral_nonlinear_diffusion_rhs(tensors, u_hat, K_hat)

    s_adv = disc.enable_nonlinear ? (disc.nonlinear_scale * disc.advection_response_scale) : 0.0
    s_diff = disc.enable_nonlinear ? (disc.nonlinear_scale * disc.diffusivity_response_scale) : 0.0

    f_tot = f_lin .- (s_adv .* f_adv) .- (s_diff .* f_diff)
    return ModalBudgetDiagnostic(f_lin, f_adv, f_diff, f_tot)
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

function _gegenbauer_weight(x::Float64, lambda::Float64)
    return max(0.0, (1.0 - x^2)^(lambda - 0.5))
end

function _reference_derivative(n::Int, lambda::Float64, x::Float64; eps::Float64 = 1e-6)
    xp = min(1.0, x + eps)
    xm = max(-1.0, x - eps)
    return (_gegenbauerC(n, lambda, xp) - _gegenbauerC(n, lambda, xm)) / (2eps)
end

function _modal_values(n_modes::Int, lambda::Float64, x::Vector{Float64})
    V = zeros(length(x), n_modes)
    for q in eachindex(x)
        for n in 1:n_modes
            V[q, n] = _gegenbauerC(n - 1, lambda, x[q])
        end
    end
    return V
end

function _modal_derivative_values(n_modes::Int, lambda::Float64, x::Vector{Float64}, H::Float64)
    dV = zeros(length(x), n_modes)
    inv_dzdx = 2.0 / H
    for q in eachindex(x)
        for n in 1:n_modes
            dV[q, n] = _reference_derivative(n - 1, lambda, x[q]) * inv_dzdx
        end
    end
    return dV
end

"""
Compute Gegenbauer triple-product tensors for nonlinear modal projections.

The tensors are defined over physical height with Gegenbauer weight:
 - triple[k,i,j] = <C_k, C_i*C_j>_lambda
 - advection[k,i,j] = <C_k, C_i*dC_j/dz>_lambda
 - diffusion_flux[k,i,j] = <-dC_k/dz, C_i*dC_j/dz>_lambda
"""
function precompute_nonlinear_tensors(disc::SpectralBLGalerkin; n_quad::Int = 256)
    n_modes = disc.n_modes
    lambda = disc.lambda
    H = disc.H

    x, w = _quad_nodes_weights(n_quad)
    V = _modal_values(n_modes, lambda, x)
    dV = _modal_derivative_values(n_modes, lambda, x, H)

    triple = zeros(n_modes, n_modes, n_modes)
    advection = zeros(n_modes, n_modes, n_modes)
    diffusion_flux = zeros(n_modes, n_modes, n_modes)

    dzdx = 0.5 * H
    for q in eachindex(x)
        xq = x[q]
        wq = _gegenbauer_weight(xq, lambda) * w[q] * dzdx
        for k in 1:n_modes
            Ck = V[q, k]
            dCk = dV[q, k]
            for i in 1:n_modes
                Ci = V[q, i]
                for j in 1:n_modes
                    Cj = V[q, j]
                    dCj = dV[q, j]

                    triple[k, i, j] += Ck * Ci * Cj * wq
                    advection[k, i, j] += Ck * Ci * dCj * wq
                    diffusion_flux[k, i, j] += -dCk * Ci * dCj * wq
                end
            end
        end
    end

    return SpectralNonlinearTensors(triple, advection, diffusion_flux)
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

    # Galerkin inner-product projection of a linear diffusive prototype operator.
    # M * da/dt = -K * a  =>  da/dt = L * a, where L = -(M \ K)
    L = _linear_modal_operator(disc)

    tensors = disc.enable_nonlinear ? precompute_nonlinear_tensors(disc) : nothing
    adv_scale, diff_scale = _closure_nonlinear_strength(closure)
    adv_scale *= disc.nonlinear_scale * disc.advection_response_scale
    diff_scale *= disc.nonlinear_scale * disc.diffusivity_response_scale

    eqs = Vector{Equation}(undef, disc.n_modes)
    for n in 1:disc.n_modes
        rhs = zero(a[1])
        for k in 1:disc.n_modes
            rhs += L[n, k] * a[k]
        end

        if disc.enable_nonlinear && !(tensors === nothing)
            for i in 1:disc.n_modes
                for j in 1:disc.n_modes
                    rhs -= adv_scale * tensors.advection[n, i, j] * a[i] * a[j]
                end
            end
            for i in 1:disc.n_modes
                for j in 1:disc.n_modes
                    rhs -= diff_scale * tensors.diffusion_flux[n, i, j] * a[i] * a[j]
                end
            end
        end

        eqs[n] = ModelingToolkit.Differential(t)(a[n]) ~ rhs
    end

    ModelingToolkit.@named modal_ode = ModelingToolkit.ODESystem(eqs, t)
    sys = structural_simplify(modal_ode)

    prob = ODEProblem(sys, ones(disc.n_modes), tspan)
    return solve(prob, solver; kwargs...)
end
