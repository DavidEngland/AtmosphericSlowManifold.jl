# src/System/SurfaceBoundary.jl
"""Return the default symbolic surface flux for a closure/state pair."""
default_surface_flux(closure::AbstractClosure, state::ManifoldState) = surface_flux(closure, state)

"""
Return the default kinematic surface sensible heat flux [K m/s].

Priority order:
1. `state.q_theta` if available.
2. `state.u_star * (state.theta - state.theta_s)` if `theta_s` exists.
3. Zero flux fallback.
"""
@inline function default_surface_heat_flux(::AbstractClosure, state::ManifoldState)
	if hasproperty(state, :q_theta)
		return getproperty(state, :q_theta)
	elseif hasproperty(state, :theta_s)
		return state.u_star * (state.theta - getproperty(state, :theta_s))
	end
	return zero(state.u_star)
end

"""
Construct default Neumann-ready surface boundary flux tuple.
"""
@inline function surface_boundary_conditions(closure::AbstractClosure, state::ManifoldState)
	tau_s = default_surface_flux(closure, state)
	q_s = default_surface_heat_flux(closure, state)
	return (momentum_flux = tau_s, heat_flux = q_s)
end
