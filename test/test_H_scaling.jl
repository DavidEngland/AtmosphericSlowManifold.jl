#!/usr/bin/env julia
# test/test_H_scaling.jl
#
# Diagnoses whether _mass_matrix (x-measure, no H) and _stiffness_matrix
# (z-measure, includes dzdx = H/2) are consistently scaled in
# SpectralBLGalerkin.jl, by calling the ACTUAL module functions directly
# (not a reimplementation) -- this only has diagnostic value if it exercises
# your real code.
#
# --- ADJUST THIS BLOCK to your real module path ----------------------------
# I'm guessing the top-level module is `AtmosphericSlowManifold` from the
# repo directory name; if SpectralBLGalerkin.jl is nested in a submodule
# (e.g. `module Discretization ... end`), extend the qualification below.
# Underscore-prefixed functions are NOT exported but ARE reachable by full
# qualification -- Julia doesn't enforce privacy on the naming convention.

using AtmosphericSlowManifold
using Test
const _mass_matrix      = AtmosphericSlowManifold._mass_matrix
const _stiffness_matrix = AtmosphericSlowManifold._stiffness_matrix
const _gegenbauerC      = AtmosphericSlowManifold._gegenbauerC
# e.g. if nested:
#   const _mass_matrix = AtmosphericSlowManifold.Discretization._mass_matrix
# -----------------------------------------------------------------------------

using LinearAlgebra
using Printf

const lambda   = 0.75
const n_modes  = 12
const Km       = 1.0    # constant diffusivity, m^2/s
const F0       = 0.05   # prescribed bottom flux magnitude

function basis_at_z(n_modes, lambda, z, H)
    x = clamp(2.0 * z / H - 1.0, -1.0, 1.0)
    return [_gegenbauerC(n - 1, lambda, x) for n in 1:n_modes]
end

reconstruct(a, n_modes, lambda, z, H) = dot(a, basis_at_z(n_modes, lambda, z, H))

# =============================================================================
# TEST 1 -- steady-state diffusion with compatible Neumann fluxes. Equal
# physical flux at both boundaries satisfies the pure-Neumann compatibility
# condition and admits a linear steady profile.
# =============================================================================
println("="^78)
println("TEST 1: steady-state Neumann diffusion slope (validates K and b)")
println("="^78)
println("theta(z) = theta(0) - (F0/Km) z  is the analytical target.")
println("Expect slope_numeric to match slope_analytic AND be identical across H.\n")

for H in (100.0, 1000.0)
    K = _stiffness_matrix(n_modes, lambda, H)     # default constant K_m profile
    phi_bottom = basis_at_z(n_modes, lambda, 0.0, H)
    phi_top = basis_at_z(n_modes, lambda, H, H)
    b = F0 .* (phi_bottom .- phi_top)

    # K has a null space (pure-Neumann constant mode) -> minimum-norm solve
    a = pinv(K) * b

    z_lo, z_hi = 0.1H, 0.9H
    slope_numeric = (reconstruct(a, n_modes, lambda, z_hi, H) -
                      reconstruct(a, n_modes, lambda, z_lo, H)) / (z_hi - z_lo)
    slope_analytic = -F0 / Km

    @printf("H = %8.1f   slope_numeric = %+.6e   slope_analytic = %+.6e   ratio = %.4f\n",
            H, slope_numeric, slope_analytic, slope_numeric / slope_analytic)
end
println()

# =============================================================================
# TEST 2 -- generalized eigenvalue check K v = lambda M v. This IS where an
# M-vs-K measure mismatch shows up: decay rates of the system are eigenvalues
# of M^{-1} K. Exact linear algebra, no time-stepping error to worry about.
# =============================================================================
println("="^78)
println("TEST 2: generalized eigenvalues of (K, M) -- validates M vs K scaling")
println("="^78)
println("Physical (Fourier-basis) analytical decay rate for the gravest non-zero")
println("mode under homogeneous Neumann BCs: lambda_1 = Km * (pi/H)^2.")
println("The Gegenbauer basis won't match this exactly (different basis functions),")
println("but the ratio numeric/analytic should be roughly CONSTANT across H if M")
println("and K share a consistent measure. If it scales by ~H or ~1/H between the")
println("two runs below, that's the missing dz/dx = H/2 factor in _mass_matrix.\n")

for H in (100.0, 1000.0)
    M = _mass_matrix(n_modes, lambda, H)
    K = _stiffness_matrix(n_modes, lambda, H)

    evals = sort(real.(eigvals(K, M)))
    scale = maximum(abs.(evals))
    nz = filter(v -> abs(v) > 1e-8 * scale, evals)
    lambda1_numeric = isempty(nz) ? NaN : nz[1]
    lambda1_analytic = Km * (pi / H)^2

    @printf("H = %8.1f   gravest_nonzero_eig = %+.6e   analytic(n=1) = %+.6e   ratio = %.4f\n",
            H, lambda1_numeric, lambda1_analytic, lambda1_numeric / lambda1_analytic)
end
println()
println("The numeric/analytic ratio should now remain constant across H.")

@testset "Spectral physical-domain metric scaling" begin
    H1, H2 = 100.0, 1000.0
    M1 = _mass_matrix(n_modes, lambda, H1)
    M2 = _mass_matrix(n_modes, lambda, H2)
    K1 = _stiffness_matrix(n_modes, lambda, H1)
    K2 = _stiffness_matrix(n_modes, lambda, H2)

    @test M2 ≈ (H2 / H1) .* M1
    @test K2 ≈ (H1 / H2) .* K1

    rates1 = sort(real.(eigvals(K1, M1)))
    rates2 = sort(real.(eigvals(K2, M2)))
    scale1 = maximum(abs, rates1)
    scale2 = maximum(abs, rates2)
    lambda1 = first(filter(value -> abs(value) > 1e-8 * scale1, rates1))
    lambda2 = first(filter(value -> abs(value) > 1e-8 * scale2, rates2))
    @test lambda2 / lambda1 ≈ (H1 / H2)^2

    boundary1 = F0 .* basis_at_z(n_modes, lambda, 0.0, H1)
    boundary2 = F0 .* basis_at_z(n_modes, lambda, 0.0, H2)
    @test boundary2 == boundary1
    @test (M2 \ boundary2) ≈ (H1 / H2) .* (M1 \ boundary1)
    @test_throws ArgumentError _mass_matrix(n_modes, lambda, 0.0)
end

@testset "Spectral nonlinear weak-load scaling" begin
    H1, H2 = 100.0, 1000.0
    modes = 6
    state = [1.0, -0.2, 0.1, 0.05, -0.02, 0.01]
    disc1 = SpectralBLGalerkin(n_modes=modes, lambda=lambda, H=H1)
    disc2 = SpectralBLGalerkin(n_modes=modes, lambda=lambda, H=H2)
    tensors1 = precompute_nonlinear_tensors(disc1; n_quad=128)
    tensors2 = precompute_nonlinear_tensors(disc2; n_quad=128)
    M1 = _mass_matrix(modes, lambda, H1; n_quad=128)
    M2 = _mass_matrix(modes, lambda, H2; n_quad=128)

    advection1 = AtmosphericSlowManifold._spectral_nonlinear_advection_rhs(tensors1, state, 1.0)
    advection2 = AtmosphericSlowManifold._spectral_nonlinear_advection_rhs(tensors2, state, 1.0)
    diffusion1 = AtmosphericSlowManifold._spectral_nonlinear_diffusion_rhs(tensors1, state, 1.0)
    diffusion2 = AtmosphericSlowManifold._spectral_nonlinear_diffusion_rhs(tensors2, state, 1.0)

    @test advection2 ≈ advection1
    @test diffusion2 ≈ (H1 / H2) .* diffusion1
    @test (M2 \ advection2) ≈ (H1 / H2) .* (M1 \ advection1)
    @test (M2 \ diffusion2) ≈ (H1 / H2)^2 .* (M1 \ diffusion1)

    positive_state = [2.0, 0.1, zeros(modes - 2)...]
    dissipative_load = AtmosphericSlowManifold._spectral_nonlinear_diffusion_rhs(
        tensors1, positive_state, 1.0
    )
    @test dot(positive_state, dissipative_load) <= 0.0
end

@testset "Spectral nonlinear weak-load Jacobian" begin
    modes = 6
    state = [1.0, -0.2, 0.1, 0.05, -0.02, 0.01]
    adv_scale, diff_scale = 0.7, 1.3
    step = 1e-6
    disc = SpectralBLGalerkin(n_modes=modes, lambda=lambda, H=100.0)
    tensors = precompute_nonlinear_tensors(disc; n_quad=128)
    jacobian = zeros(modes, modes)
    AtmosphericSlowManifold._spectral_nonlinear_jacobian!(
        jacobian, tensors, state, adv_scale, diff_scale
    )

    load(candidate) = AtmosphericSlowManifold._spectral_nonlinear_rhs(
        tensors, candidate, adv_scale, diff_scale
    )
    finite_difference = zeros(modes, modes)
    for column in 1:modes
        perturbation = zeros(modes)
        perturbation[column] = step
        finite_difference[:, column] .= (
            load(state .+ perturbation) .- load(state .- perturbation)
        ) ./ (2step)
    end

    @test jacobian ≈ finite_difference rtol=1e-8 atol=1e-9
end