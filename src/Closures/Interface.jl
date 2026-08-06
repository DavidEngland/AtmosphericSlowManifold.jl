abstract type AbstractClosure end

function eddy_momentum(::AbstractClosure, ::ManifoldState)
    throw(MethodError(eddy_momentum, (AbstractClosure, ManifoldState)))
end

function eddy_heat(::AbstractClosure, ::ManifoldState)
    throw(MethodError(eddy_heat, (AbstractClosure, ManifoldState)))
end

function surface_flux(::AbstractClosure, ::ManifoldState)
    throw(MethodError(surface_flux, (AbstractClosure, ManifoldState)))
end
