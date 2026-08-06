"""Return the default symbolic surface flux for a closure/state pair."""
default_surface_flux(closure::AbstractClosure, state::ManifoldState) = surface_flux(closure, state)
