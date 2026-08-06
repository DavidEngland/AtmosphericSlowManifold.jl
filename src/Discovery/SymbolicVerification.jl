function _try_eval(expr::Num, substitutions::Dict{Num, <:Real})
    value = Symbolics.value(Symbolics.substitute(expr, substitutions))
    return value isa Number ? float(value) : NaN
end

"""Run basic positivity and fold-readiness checks for a closure."""
function verify_closure(
    closure::AbstractClosure;
    state::ManifoldState = ManifoldState(),
    substitutions::Dict{Num, <:Real} = Dict{Num, Float64}()
)
    km = eddy_momentum(closure, state)
    kh = eddy_heat(closure, state)
    flux = surface_flux(closure, state)

    km_val = isempty(substitutions) ? NaN : _try_eval(km, substitutions)
    kh_val = isempty(substitutions) ? NaN : _try_eval(kh, substitutions)

    return Dict(
        :km_expr => km,
        :kh_expr => kh,
        :flux_expr => flux,
        :km_positive => isnan(km_val) ? missing : km_val > 0,
        :kh_positive => isnan(kh_val) ? missing : kh_val > 0,
    )
end
