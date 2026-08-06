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
    z0::Num
end

function ManifoldState(; name::Symbol = :manifold)
    ModelingToolkit.@variables eta1 eta2 eta3 r omega chi pi_g lambdamin u v theta q u_star z0
    return ManifoldState(eta1, eta2, eta3, r, omega, chi, pi_g, lambdamin, u, v, theta, q, u_star, z0)
end
