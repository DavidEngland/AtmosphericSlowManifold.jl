struct WSINDyClosure <: AbstractClosure
    km_expr::Num
    kh_expr::Num
    flux_expr::Num
end

eddy_momentum(c::WSINDyClosure, ::ManifoldState) = c.km_expr
eddy_heat(c::WSINDyClosure, ::ManifoldState) = c.kh_expr
surface_flux(c::WSINDyClosure, ::ManifoldState) = c.flux_expr
