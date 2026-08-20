# src/Manifold/ManifoldState.jl
struct ManifoldState
    eta1::Num
    eta2::Num
    eta3::Num
    r::Num
    omega::Num
    chi::Num
    pi_g::Num
    lambdamin::Num
    u::Num
    v::Num
    theta::Num
    q::Num
    u_star::Num
    z::Num
    z0::Num
end

function ManifoldState(;
    name::Symbol = :manifold,
    eta1 = nothing, eta2 = nothing, eta3 = nothing, r = nothing, omega = nothing,
    chi = nothing, pi_g = nothing, lambdamin = nothing, u = nothing, v = nothing,
    theta = nothing, q = nothing, u_star = nothing, z = nothing, z0 = nothing,
)
    ModelingToolkit.@variables sym_eta1 sym_eta2 sym_eta3 sym_r sym_omega sym_chi sym_pi_g sym_lambdamin sym_u sym_v sym_theta sym_q sym_u_star sym_z sym_z0

    # Numeric overrides fall back to fresh symbolic variables when omitted.
    pick(val, sym) = val === nothing ? sym : Num(val)

    return ManifoldState(
        pick(eta1, sym_eta1), pick(eta2, sym_eta2), pick(eta3, sym_eta3),
        pick(r, sym_r), pick(omega, sym_omega), pick(chi, sym_chi),
        pick(pi_g, sym_pi_g), pick(lambdamin, sym_lambdamin),
        pick(u, sym_u), pick(v, sym_v), pick(theta, sym_theta),
        pick(q, sym_q), pick(u_star, sym_u_star), pick(z, sym_z), pick(z0, sym_z0),
    )
end
