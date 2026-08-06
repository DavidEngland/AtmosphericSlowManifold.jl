struct MOSTClosure <: AbstractClosure
    kappa::Float64
    pr_t::Float64
    cm_stable::Float64
end

function MOSTClosure(; kappa::Float64 = 0.41, pr_t::Float64 = 0.74, cm_stable::Float64 = 4.7)
    return MOSTClosure(kappa, pr_t, cm_stable)
end

function eddy_momentum(c::MOSTClosure, m::ManifoldState)
    # Minimal symbolic baseline for neutral-to-stable behavior.
    return (c.kappa * m.u_star * m.z0) / (1 + c.cm_stable * m.r)
end

eddy_heat(c::MOSTClosure, m::ManifoldState) = eddy_momentum(c, m) / c.pr_t
surface_flux(::MOSTClosure, m::ManifoldState) = m.u_star^2
