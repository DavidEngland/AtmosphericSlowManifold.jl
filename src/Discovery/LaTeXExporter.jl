# s
using Symbolics

export to_latex, latex_term_table, latex_site_summary_table

const _GREEK_MAP = Dict(
    "alpha" => "\\alpha",
    "beta" => "\\beta",
    "gamma" => "\\gamma",
    "delta" => "\\delta",
    "epsilon" => "\\epsilon",
    "theta" => "\\theta",
    "kappa" => "\\kappa",
    "lambda" => "\\lambda",
    "mu" => "\\mu",
    "nu" => "\\nu",
    "xi" => "\\xi",
    "pi" => "\\pi",
    "rho" => "\\rho",
    "sigma" => "\\sigma",
    "tau" => "\\tau",
    "phi" => "\\phi",
    "psi" => "\\psi",
    "omega" => "\\omega",
)

function _num_to_latex(x::Real)
    if iszero(x)
        return "0"
    end

    if isinteger(x)
        return string(Int(round(x)))
    end

    xs = string(round(Float64(x); sigdigits = 6))
    if occursin("e", xs)
        parts = split(xs, "e")
        base = parts[1]
        expo = parse(Int, parts[2])
        return string(base, "\\times 10^{", expo, "}")
    end
    return xs
end

function _latex_symbol(sym::Symbol)
    s = String(sym)
    if haskey(_GREEK_MAP, lowercase(s))
        return _GREEK_MAP[lowercase(s)]
    elseif length(s) == 1
        return s
    end

    escaped = replace(s, "_" => "\\_")
    return string("\\mathrm{", escaped, "}")
end

function _feature_to_latex(f::StateVariable)
    return _latex_symbol(f.name)
end

function _feature_to_latex(f::DiagnosticVariable)
    return _latex_symbol(f.name)
end

function _feature_to_latex(f::SpatialDerivative)
    var = _latex_symbol(f.variable)
    z = "z"
    if f.order <= 1
        return string("\\frac{\\partial ", var, "}{\\partial ", z, "}")
    end
    return string("\\frac{\\partial^{", f.order, "} ", var, "}{\\partial ", z, "^{", f.order, "}}")
end

function _basis_to_latex(b::BasisOperator)
    base = _feature_to_latex(b.feature)
    p = b.power
    if abs(p - 1.0) <= 1e-12
        return base
    end

    ptxt = isinteger(p) ? string(Int(round(p))) : _num_to_latex(p)
    return string("\\left(", base, "\\right)^{", ptxt, "}")
end

function _term_magnitude_to_latex(term::OperatorTerm)
    c = abs(Float64(term.coefficient))
    factors = [_basis_to_latex(b) for b in term.basis]

    if isempty(factors)
        return _num_to_latex(c)
    end

    if abs(c - 1.0) <= 1e-12
        return join(factors, " ")
    end

    return string(_num_to_latex(c), " ", join(factors, " "))
end

function _rhs_to_latex(model::DiscoveredModel)
    if isempty(model.terms)
        return "0"
    end

    chunks = String[]
    for (i, term) in enumerate(model.terms)
        sign = Float64(term.coefficient) < 0 ? "-" : "+"
        mag = _term_magnitude_to_latex(term)

        if i == 1
            if sign == "-"
                push!(chunks, string("-", mag))
            else
                push!(chunks, mag)
            end
        else
            push!(chunks, string(" ", sign, " ", mag))
        end
    end
    return join(chunks)
end

"""
    to_latex(model::DiscoveredModel; independent=:t)

Create publication-ready LaTeX for discovered prognostic equations.
"""
function to_latex(model::DiscoveredModel; independent::Symbol = :t)
    lhs = string(
        "\\frac{\\partial ",
        _latex_symbol(model.target_variable),
        "}{\\partial ",
        _latex_symbol(independent),
        "}",
    )
    rhs = _rhs_to_latex(model)
    return string(lhs, " = ", rhs)
end

function _symbolics_text_to_latex(expr)
    s = string(expr)
    s = replace(s, "θ" => "\\theta")
    s = replace(s, "ν" => "\\nu")
    s = replace(s, "ξ" => "\\xi")
    s = replace(s, "κ" => "\\kappa")
    return replace(s, "_" => "\\_")
end

"""
    to_latex(closure::WSINDyClosure; field=:momentum)

Return LaTeX for one closure expression (`:momentum`, `:heat`, or `:flux`).
"""
function to_latex(closure::WSINDyClosure; field::Symbol = :momentum)
    expr = if field == :momentum
        closure.km_expr
    elseif field == :heat
        closure.kh_expr
    elseif field == :flux
        closure.flux_expr
    else
        throw(ArgumentError("Unknown WSINDyClosure field $(field). Use :momentum, :heat, or :flux."))
    end
    return _symbolics_text_to_latex(expr)
end

"""
    latex_term_table(model::DiscoveredModel; r2=nothing, residual_norm=model.residual_norm)

Generate a LaTeX tabular environment for discovered terms and coefficients.
"""
function latex_term_table(
    model::DiscoveredModel;
    r2::Union{Nothing, Real} = nothing,
    residual_norm::Real = model.residual_norm,
)
    lines = String[]
    push!(lines, "\\begin{tabular}{r l r}")
    push!(lines, "\\hline")
    push!(lines, "Term & Basis & Coefficient \\\\")
    push!(lines, "\\hline")

    for (i, term) in enumerate(model.terms)
        basis = isempty(term.basis) ? "1" : join([_basis_to_latex(b) for b in term.basis], " ")
        coeff = _num_to_latex(Float64(term.coefficient))
        push!(lines, string(i, " & \$", basis, "\$ & \$", coeff, "\$ \\\\"))
    end

    push!(lines, "\\hline")
    if r2 !== nothing
        push!(lines, string("\\multicolumn{3}{l}{\$R^{2}=", _num_to_latex(Float64(r2)), "\$} \\\\"))
    end
    push!(lines, string("\\multicolumn{3}{l}{Residual norm=", _num_to_latex(Float64(residual_norm)), "} \\\\"))
    push!(lines, "\\hline")
    push!(lines, "\\end{tabular}")
    return join(lines, "\n")
end

"""
    latex_site_summary_table(models; r2_by_site=Dict(), residual_by_site=Dict())

Create a publication table comparing fit quality across sites.
"""
function latex_site_summary_table(
    models::AbstractDict{Symbol, <:DiscoveredModel};
    r2_by_site::Dict{Symbol, <:Real} = Dict{Symbol, Float64}(),
    residual_by_site::Dict{Symbol, <:Real} = Dict{Symbol, Float64}(),
)
    lines = String[]
    push!(lines, "\\begin{tabular}{l r r r}")
    push!(lines, "\\hline")
    push!(lines, "Site & Active terms & \$R^{2}\$ & Residual norm \\\\")
    push!(lines, "\\hline")

    for site in sort(collect(keys(models)); by = String)
        m = models[site]
        n_active = m.sparsity_level
        r2 = haskey(r2_by_site, site) ? _num_to_latex(Float64(r2_by_site[site])) : "-"
        rn = haskey(residual_by_site, site) ? _num_to_latex(Float64(residual_by_site[site])) : _num_to_latex(m.residual_norm)
        push!(lines, string(_latex_symbol(site), " & ", n_active, " & ", r2, " & ", rn, " \\\\"))
    end

    push!(lines, "\\hline")
    push!(lines, "\\end{tabular}")
    return join(lines, "\n")
end