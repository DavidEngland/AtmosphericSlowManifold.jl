abstract type AbstractClosure end

# Compatibility alias used by external notes/specs.
const AbstractAtmosphericClosure = AbstractClosure

function eddy_momentum(::AbstractClosure, ::ManifoldState)
    throw(MethodError(eddy_momentum, (AbstractClosure, ManifoldState)))
end

function eddy_heat(::AbstractClosure, ::ManifoldState)
    throw(MethodError(eddy_heat, (AbstractClosure, ManifoldState)))
end

function surface_flux(::AbstractClosure, ::ManifoldState)
    throw(MethodError(surface_flux, (AbstractClosure, ManifoldState)))
end

function evaluate_diffusivity_profile!(::AbstractVector, ::AbstractClosure, ::AbstractVector)
    throw(MethodError(evaluate_diffusivity_profile!, (AbstractVector, AbstractClosure, AbstractVector)))
end

function evaluate_heat_diffusivity_profile!(::AbstractVector, ::AbstractClosure, ::AbstractVector)
    throw(MethodError(evaluate_heat_diffusivity_profile!, (AbstractVector, AbstractClosure, AbstractVector)))
end
