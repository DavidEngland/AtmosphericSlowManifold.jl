# src/Discretization/Backends/SpectralBLGalerkin.jl
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
    triple::Array{Float64,3}
    advection::Array{Float64,3}
    diffusion_flux::Array{Float64,3}
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

"""
Pre-allocated workspace for closure-driven diffusivity coupling in spectral space.
"""
struct BoundaryLayerWorkspace{T<:AbstractFloat,C<:AbstractClosure}
    closure::C
    z_grid::Vector{T}
    K_m_buffer::Vector{T}
    K_h_buffer::Vector{T}
    dK_dz_buffer::Vector{T}
end

function build_boundary_layer_workspace(disc::SpectralBLGalerkin, closure::AbstractClosure; z_min::Float64=2.0)
    z_grid = collect(range(z_min, disc.H; length=disc.n_modes))
    T = Float64
    return BoundaryLayerWorkspace{T,typeof(closure)}(
        closure,
        T.(z_grid),
        zeros(T, disc.n_modes),
        zeros(T, disc.n_modes),
        zeros(T, disc.n_modes),
    )
end

function _fill_gradient!(dK_dz::Vector{Float64}, K::Vector{Float64}, z::Vector{Float64})
    n = length(K)
    (n == length(dK_dz) && n == length(z)) || throw(DimensionMismatch("Gradient buffers must have matching lengths."))
    if n == 0
        return dK_dz
    elseif n == 1
        dK_dz[1] = 0.0
        return dK_dz
    end

    @inbounds begin
        dK_dz[1] = (K[2] - K[1]) / (z[2] - z[1])
        for i in 2:(n-1)
            dK_dz[i] = (K[i+1] - K[i-1]) / (z[i+1] - z[i-1])
        end
        dK_dz[n] = (K[n] - K[n-1]) / (z[n] - z[n-1])
    end
    return dK_dz
end

function update_diffusivity_buffers!(workspace::BoundaryLayerWorkspace{Float64})
    closure = workspace.closure
    if closure isa PhysicalSimilarityClosure
        evaluate_diffusivity_profile!(workspace.K_m_buffer, closure, workspace.z_grid)
        evaluate_heat_diffusivity_profile!(workspace.K_h_buffer, closure, workspace.z_grid)
    else
        adv_scale, diff_scale = _closure_nonlinear_strength(closure)
        fill!(workspace.K_m_buffer, abs(adv_scale))
        fill!(workspace.K_h_buffer, abs(diff_scale))
    end
    _fill_gradient!(workspace.dK_dz_buffer, workspace.K_m_buffer, workspace.z_grid)
    return workspace
end

@inline function _mean_abs(v::Vector{Float64})
    n = length(v)
    n == 0 && return 0.0
    acc = 0.0
    @inbounds for i in eachindex(v)
        acc += abs(v[i])
    end
    return acc / n
end

@inline function _interp1d(x_target::Float64, x_grid::Vector{Float64}, y_grid::Vector{Float64})
    n = length(x_grid)
    n == 0 && return 0.0
    n == 1 && return y_grid[1]
    if x_target <= x_grid[1]
        return y_grid[1]
    elseif x_target >= x_grid[end]
        return y_grid[end]
    end
    idx = searchsortedlast(x_grid, x_target)
    idx = clamp(idx, 1, n - 1)
    x0, x1 = x_grid[idx], x_grid[idx+1]
    y0, y1 = y_grid[idx], y_grid[idx+1]
    t = (x_target - x0) / (x1 - x0)
    return y0 + t * (y1 - y0)
end

function spectral_rhs!(
    du::Vector{Float64},
    u::Vector{Float64},
    L::Matrix{Float64},
    tensors::Union{Nothing,SpectralNonlinearTensors},
    workspace::BoundaryLayerWorkspace{Float64},
    disc::SpectralBLGalerkin,
    adv_scale::Float64,
    diff_scale::Float64,
)
    mul!(du, L, u)
    if !disc.enable_nonlinear || tensors === nothing
        return du
    end

    n = disc.n_modes
    @inbounds for k in 1:n
        adv = 0.0
        diff = 0.0
        for i in 1:n
            ui = u[i]
            for j in 1:n
                uij = ui * u[j]
                adv += tensors.advection[k, i, j] * uij
                diff += tensors.diffusion_flux[k, i, j] * uij
            end
        end
        du[k] -= adv_scale * adv + diff_scale * diff
    end
    return du
end

function SpectralBLGalerkin(
    ;
    n_modes::Int=12,
    lambda::Float64=0.75,
    H::Float64=3000.0,
    enable_nonlinear::Bool=true,
    nonlinear_scale::Float64=1.0,
    advection_response_scale::Float64=1.0,
    diffusivity_response_scale::Float64=1.0,
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

function _linear_modal_operator(
    disc::SpectralBLGalerkin,
    K_m_profile::Union{Nothing,Vector{Float64}}=nothing,
    z_grid::Union{Nothing,Vector{Float64}}=nothing,
)
    M = _mass_matrix(disc.n_modes, disc.lambda)
    K = _stiffness_matrix(disc.n_modes, disc.lambda, disc.H, K_m_profile, z_grid)
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
    z_grid::Union{Nothing,Vector{Float64}}=nothing,
)
    length(u_hat) == disc.n_modes || throw(ArgumentError("u_hat length must match disc.n_modes."))

    Lop = _linear_modal_operator(disc, K_hat, z_grid)
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

function _modal_basis_at_z(n_modes::Int, lambda::Float64, z::Float64, H::Float64)
    0.0 <= z <= H || throw(ArgumentError("Basis evaluation height must lie in [0, H]."))
    x = 2.0 * z / H - 1.0
    return [_gegenbauerC(n - 1, lambda, x) for n in 1:n_modes]
end

function _add_surface_heat_forcing!(
    du::AbstractVector,
    u::AbstractVector,
    forcing::SurfaceForcing,
    phi_z1::Vector{Float64},
    transfer_velocity::Float64,
    t::Real,
)
    surface_temperature = interp_forcing(forcing, :Ts, t)
    reference_temperature = dot(phi_z1, u)
    heat_flux = transfer_velocity * (surface_temperature - reference_temperature)
    @inbounds for mode in eachindex(phi_z1)
        du[mode] += phi_z1[mode] * heat_flux
    end
    return du
end

function _add_surface_heat_jacobian!(
    jacobian::AbstractMatrix,
    phi_z1::Vector{Float64},
    transfer_velocity::Float64,
)
    @inbounds for row in eachindex(phi_z1), column in eachindex(phi_z1)
        jacobian[row, column] -= transfer_velocity * phi_z1[row] * phi_z1[column]
    end
    return jacobian
end

function _x_to_z(x::Float64, H::Float64)
    return 0.5 * H * (x + 1.0)
end

function _quad_nodes_weights(n::Int)
    xs = range(-1.0, 1.0; length=n)
    dx = step(xs)
    ws = fill(dx, n)
    ws[1] *= 0.5
    ws[end] *= 0.5
    return collect(xs), ws
end

function _gegenbauer_weight(x::Float64, lambda::Float64)
    return max(0.0, (1.0 - x^2)^(lambda - 0.5))
end

function _reference_derivative(n::Int, lambda::Float64, x::Float64; eps::Float64=1e-6)
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
function precompute_nonlinear_tensors(disc::SpectralBLGalerkin; n_quad::Int=256)
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

function _mass_matrix(n_modes::Int, lambda::Float64; n_quad::Int=256)
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

function _stiffness_matrix(
    n_modes::Int,
    lambda::Float64,
    H::Float64,
    K_m_profile::Union{Nothing,Vector{Float64}}=nothing,
    z_grid::Union{Nothing,Vector{Float64}}=nothing;
    n_quad::Int=256,
)
    x, w = _quad_nodes_weights(n_quad)
    K = zeros(n_modes, n_modes)
    dzdx = 0.5 * H
    inv_dzdx = 1.0 / dzdx
    eps = 1e-6

    Km_quad = ones(length(x))
    if K_m_profile !== nothing && z_grid !== nothing && !isempty(K_m_profile)
        for q in eachindex(x)
            zq = _x_to_z(x[q], H)
            Km_quad[q] = max(0.0, _interp1d(zq, z_grid, K_m_profile))
        end
    end

    for i in 1:n_modes
        for j in 1:n_modes
            acc = 0.0
            ni = i - 1
            nj = j - 1
            for q in eachindex(x)
                xq = x[q]
                weight = (1.0 - xq^2)^(lambda - 0.5)

                dCi_dx = (_gegenbauerC(ni, lambda, min(1.0, xq + eps)) - _gegenbauerC(ni, lambda, max(-1.0, xq - eps))) / (2 * eps)
                dCj_dx = (_gegenbauerC(nj, lambda, min(1.0, xq + eps)) - _gegenbauerC(nj, lambda, max(-1.0, xq - eps))) / (2 * eps)
                dCi_dz = dCi_dx * inv_dzdx
                dCj_dz = dCj_dx * inv_dzdx

                acc += Km_quad[q] * dCi_dz * dCj_dz * weight * dzdx * w[q]
            end
            K[i, j] = acc
        end
    end
    return K
end

"""
Solve the single-field spectral system.

When `surface_forcing` is provided, the modal state is interpreted as potential
temperature in Kelvin. `bulk_transfer_coeff` and `reference_wind_speed` define
the bulk transfer velocity, and `z1` is the reference measurement height.
"""
function dispatch_solve(
    disc::SpectralBLGalerkin,
    pde_sys::PDESystem,
    closure::AbstractClosure,
    tspan::Tuple{Float64,Float64};
    solver=Rodas5P(autodiff=true),
    u0=nothing,
    reltol::Real=1e-6,
    abstol::Real=1e-8,
    maxiters::Int=1_000_000,
    surface_forcing::Union{Nothing,SurfaceForcing}=nothing,
    bulk_transfer_coeff::Float64=1.5e-3,
    reference_wind_speed::Float64=5.0,
    z1::Union{Nothing,Float64}=nothing,
    kwargs...
)
    workspace = build_boundary_layer_workspace(disc, closure)
    update_diffusivity_buffers!(workspace)

    bulk_transfer_coeff >= 0.0 || throw(ArgumentError("bulk_transfer_coeff must be nonnegative."))
    reference_wind_speed >= 0.0 || throw(ArgumentError("reference_wind_speed must be nonnegative."))
    z1_actual = z1 === nothing ? first(workspace.z_grid) : z1
    phi_z1 = _modal_basis_at_z(disc.n_modes, disc.lambda, z1_actual, disc.H)
    transfer_velocity = bulk_transfer_coeff * reference_wind_speed

    M = _mass_matrix(disc.n_modes, disc.lambda)
    K = _stiffness_matrix(disc.n_modes, disc.lambda, disc.H, workspace.K_m_buffer, workspace.z_grid)
    tensors = disc.enable_nonlinear ? precompute_nonlinear_tensors(disc) : nothing
    adv_scale = _mean_abs(workspace.K_m_buffer)
    diff_scale = _mean_abs(workspace.K_h_buffer)
    adv_scale *= disc.nonlinear_scale * disc.advection_response_scale
    diff_scale *= disc.nonlinear_scale * disc.diffusivity_response_scale

    init = if u0 === nothing
        ones(disc.n_modes)
    else
        length(u0) == disc.n_modes || throw(DimensionMismatch("u0 length must match disc.n_modes."))
        Float64.(collect(u0))
    end

    p = (
        M=M,
        K=K,
        tensors=tensors,
        adv_scale=adv_scale,
        diff_scale=diff_scale,
        n=disc.n_modes,
        nl=zeros(Float64, disc.n_modes),
        dnl=zeros(Float64, disc.n_modes, disc.n_modes),
        M_dnl=zeros(Float64, disc.n_modes, disc.n_modes),
        surface_forcing=surface_forcing,
        phi_z1=phi_z1,
        transfer_velocity=transfer_velocity,
    )

    function rhs_mass!(du, u, p, t)
        fill!(du, 0.0)
        mul!(du, p.K, u)
        @inbounds for i in 1:p.n
            du[i] = -du[i]
        end

        if p.tensors !== nothing
            fill!(p.nl, 0.0)
            @inbounds for k in 1:p.n
                adv = 0.0
                diff = 0.0
                for i in 1:p.n
                    ui = u[i]
                    for j in 1:p.n
                        uij = ui * u[j]
                        adv += p.tensors.advection[k, i, j] * uij
                        diff += p.tensors.diffusion_flux[k, i, j] * uij
                    end
                end
                p.nl[k] = -(p.adv_scale * adv + p.diff_scale * diff)
            end

            mul!(p.nl, p.M, p.nl)
            @inbounds for i in 1:p.n
                du[i] += p.nl[i]
            end
        end

        p.surface_forcing === nothing || _add_surface_heat_forcing!(
            du,
            u,
            p.surface_forcing,
            p.phi_z1,
            p.transfer_velocity,
            t,
        )
        return du
    end

    function jac_mass!(J, u, p, t)
        fill!(J, 0.0)
        @inbounds for i in 1:p.n
            for j in 1:p.n
                J[i, j] = -p.K[i, j]
            end
        end

        if p.tensors !== nothing
            fill!(p.dnl, 0.0)
            @inbounds for k in 1:p.n
                for l in 1:p.n
                    d_adv = 0.0
                    d_diff = 0.0
                    for j in 1:p.n
                        d_adv += p.tensors.advection[k, l, j] * u[j]
                        d_adv += p.tensors.advection[k, j, l] * u[j]
                        d_diff += p.tensors.diffusion_flux[k, l, j] * u[j]
                        d_diff += p.tensors.diffusion_flux[k, j, l] * u[j]
                    end
                    p.dnl[k, l] = -(p.adv_scale * d_adv + p.diff_scale * d_diff)
                end
            end

            mul!(p.M_dnl, p.M, p.dnl)
            @inbounds for i in 1:p.n
                for j in 1:p.n
                    J[i, j] += p.M_dnl[i, j]
                end
            end
        end

        p.surface_forcing === nothing || _add_surface_heat_jacobian!(
            J,
            p.phi_z1,
            p.transfer_velocity,
        )
        return J
    end

    odef = ODEFunction(rhs_mass!; jac=jac_mass!, mass_matrix=M)
    prob = ODEProblem(odef, init, tspan, p)
    return solve(prob, solver; reltol=reltol, abstol=abstol, maxiters=maxiters, kwargs...)
end