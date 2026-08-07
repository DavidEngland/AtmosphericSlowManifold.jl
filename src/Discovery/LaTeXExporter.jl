# src/Discovery/LaTeXExporter.jl
using Latexify
using Printf
using Symbolics

export to_latex, latex_term_table, latex_site_summary_table, write_latex

const _GREEK_MAP = Dict(
    "alpha" => "\\alpha",
    "beta" => "\\beta",
    "gamma" => "\\gamma",
    "delta" => "\\delta",
    "epsilon" => "\\epsilon",
    "zeta" => "\\zeta",
    "eta" => "\\eta",
    "theta" => "\\theta",
    "iota" => "\\iota",
    "kappa" => "\\kappa",
    "lambda" => "\\lambda",
    "mu" => "\\mu",
    "nu" => "\\nu",
    "xi" => "\\xi",
    "pi" => "\\pi",
    "rho" => "\\rho",
    "sigma" => "\\sigma",
    "tau" => "\\tau",
    "upsilon" => "\\upsilon",
    "phi" => "\\phi",
    "chi" => "\\chi",
    "psi" => "\\psi",
    "omega" => "\\omega",
    "Gamma" => "\\Gamma",
    "Delta" => "\\Delta",
    "Theta" => "\\Theta",
    "Lambda" => "\\Lambda",
    "Xi" => "\\Xi",
    "Pi" => "\\Pi",
    "Sigma" => "\\Sigma",
    "Phi" => "\\Phi",
    "Psi" => "\\Psi",
    "Omega" => "\\Omega",
)

function _num_to_latex(x::Real; sigdigits::Int = 6)
    if iszero(x)
        return "0"
    elseif isinteger(x)
        return string(Int(round(x)))
    end

    formatted = @sprintf("%.*g", sigdigits, Float64(x))
    if occursin("e", formatted)
        base, expo = split(formatted, "e")
        return string(base, "\\times 10^{", parse(Int, expo), "}")
    end
    return formatted
end

function _format_symbol_part(part::AbstractString)
    if haskey(_GREEK_MAP, part)
        return _GREEK_MAP[part]
    elseif haskey(_GREEK_MAP, lowercase(part))
        return _GREEK_MAP[lowercase(part)]
    elseif isempty(part)
        return part
    elseif length(part) == 1 && isletter(first(part))
        return part
    end

    escaped = replace(part, "-" => "\\text{-}", "_" => "\\_")
    return string("\\mathrm{", escaped, "}")
end

function _latex_symbol(sym::Symbol)
    parts = split(String(sym), "_")
    head = _format_symbol_part(first(parts))
    if length(parts) == 1
        return head
    end

    subparts = join((_format_symbol_part(p) for p in parts[2:end]), "\\,")
    return string(head, "_{", subparts, "}")
end

function _feature_to_latex(f::StateVariable)
    return _latex_symbol(f.name)
end

function _feature_to_latex(f::DiagnosticVariable)
    return _latex_symbol(f.name)
end

function _feature_to_latex(f::SpatialDerivative)
    var = _latex_symbol(f.variable)
    coordinate = hasproperty(f, :coordinate) ? getfield(f, :coordinate) : :z
    coord = _latex_symbol(coordinate)
    if f.order <= 1
        return string("\\frac{\\partial ", var, "}{\\partial ", coord, "}")
    end
    return string("\\frac{\\partial^{", f.order, "} ", var, "}{\\partial ", coord, "^{", f.order, "}}")
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

function _term_to_latex(term::OperatorTerm)
    coeff = Float64(term.coefficient)
    magnitude = _num_to_latex(abs(coeff))
    factors = isempty(term.basis) ? "" : join((_basis_to_latex(b) for b in term.basis), "\\,")
    if isempty(factors)
        return magnitude
    end
    if abs(abs(coeff) - 1.0) <= 1e-12
        return factors
    end
    return string(magnitude, "\\,", factors)
end

function _rhs_to_latex(model::DiscoveredModel)
    isempty(model.terms) && return "0"

    chunks = String[]
    for (i, term) in enumerate(model.terms)
        term_str = _term_to_latex(term)
        if i == 1
            if Float64(term.coefficient) < 0
                push!(chunks, string("-", term_str))
            else
                push!(chunks, term_str)
            end
        else
            sign = Float64(term.coefficient) < 0 ? "-" : "+"
            push!(chunks, string(" ", sign, " ", term_str))
        end
    end
    return join(chunks)
end

function _wrap_equation(expr::String, environment::Symbol)
    if environment == :none
        return expr
    elseif environment in (:equation, :eq)
        return string("\\begin{equation}\n", expr, "\n\\end{equation}")
    elseif environment in (:equation_star, :eq_star, Symbol("equation*"))
        return string("\\begin{equation*}\n", expr, "\n\\end{equation*}")
    elseif environment == :align
        return string("\\begin{align}\n", expr, "\n\\end{align}")
    elseif environment in (:align_star, Symbol("align*"))
        return string("\\begin{align*}\n", expr, "\n\\end{align*}")
    else
        throw(ArgumentError("Unsupported environment $(environment)."))
    end
end

"""
    to_latex(model::DiscoveredModel; independent=:t, environment=:none)

Render a discovered model as a publication-ready LaTeX equation.
"""
function to_latex(model::DiscoveredModel; independent::Symbol = :t, environment::Symbol = :none)
    lhs = string("\\frac{\\partial ", _latex_symbol(model.target_variable), "}{\\partial ", _latex_symbol(independent), "}")
    expr = string(lhs, " = ", _rhs_to_latex(model))
    return _wrap_equation(expr, environment)
end

"""
    to_latex(models::AbstractVector{<:DiscoveredModel}; independent=:t, environment=:align)

Render multiple discovered models as an aligned LaTeX system.
"""
function to_latex(models::AbstractVector{<:DiscoveredModel}; independent::Symbol = :t, environment::Symbol = :align)
    lines = String[]
    for m in models
        lhs = string("\\frac{\\partial ", _latex_symbol(m.target_variable), "}{\\partial ", _latex_symbol(independent), "}")
        push!(lines, string(lhs, " &= ", _rhs_to_latex(m)))
    end
    body = join(lines, " " * "\\\\" * "\n")
    if environment == :none
        return body
    elseif environment == :align
        return string("\\begin{align}\n", body, "\n\\end{align}")
    elseif environment in (:align_star, Symbol("align*"))
        return string("\\begin{align*}\n", body, "\n\\end{align*}")
    else
        return body
    end
end

"""
    to_latex(closure::WSINDyClosure; field=:momentum)

Render a closure expression using Latexify for robust symbolic formatting.
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
    return String(latexify(Symbolics.simplify(expr)))
end

"""
    write_latex(path::AbstractString, content::AbstractString)

Write a LaTeX string directly to a file path, ensuring parent directories exist.
"""
function write_latex(path::AbstractString, content::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, content)
    end
    return path
end

"""
    latex_term_table(model::DiscoveredModel; r2=nothing, residual_norm=model.residual_norm, booktabs=true, longtable=false)

Generate a LaTeX `tabular` or `longtable` environment for discovered terms and coefficients.
"""
function latex_term_table(
    model::DiscoveredModel;
    r2::Union{Nothing, Real} = nothing,
    residual_norm::Real = model.residual_norm,
    booktabs::Bool = true,
    longtable::Bool = false,
)
    table_env = longtable ? "longtable" : "tabular"
    top_rule = booktabs ? "\\toprule" : "\\hline"
    mid_rule = booktabs ? "\\midrule" : "\\hline"
    bot_rule = booktabs ? "\\bottomrule" : "\\hline"

    lines = String[]
    push!(lines, string("\\begin{", table_env, "}{r l r}"))
    push!(lines, top_rule)
    push!(lines, string("Term & Basis & Coefficient ", "\\\\"))
    push!(lines, mid_rule)

    for (i, term) in enumerate(model.terms)
        basis = isempty(term.basis) ? "1" : join((_basis_to_latex(b) for b in term.basis), "\\,")
        coeff = _num_to_latex(Float64(term.coefficient))
        push!(lines, string(i, " & ", '\$', basis, '\$', " & ", '\$', coeff, '\$', " ", "\\\\"))
    end

    push!(lines, mid_rule)
    if r2 !== nothing
        push!(lines, string("\\multicolumn{3}{l}{", '\$', "R^{2} = ", _num_to_latex(Float64(r2)), '\$', "} ", "\\\\"))
    end
    push!(lines, string("\\multicolumn{3}{l}{Residual norm = ", _num_to_latex(Float64(residual_norm)), "} ", "\\\\"))
    push!(lines, bot_rule)
    push!(lines, string("\\end{", table_env, "}"))
    return join(lines, "\n")
end

"""
    latex_site_summary_table(models; r2_by_site=Dict(), residual_by_site=Dict(), booktabs=true)

Create a publication table comparing fit quality across sites.
"""
function latex_site_summary_table(
    models::AbstractDict{Symbol, <:DiscoveredModel};
    r2_by_site::Dict{Symbol, <:Real} = Dict{Symbol, Float64}(),
    residual_by_site::Dict{Symbol, <:Real} = Dict{Symbol, Float64}(),
    booktabs::Bool = true,
)
    top_rule = booktabs ? "\\toprule" : "\\hline"
    mid_rule = booktabs ? "\\midrule" : "\\hline"
    bot_rule = booktabs ? "\\bottomrule" : "\\hline"

    lines = String[]
    push!(lines, "\\begin{tabular}{l r r r}")
    push!(lines, top_rule)
    push!(lines, string("Site & Active terms & ", '\$', "R^{2}", '\$', " & Residual norm ", "\\\\"))
    push!(lines, mid_rule)

    for site in sort(collect(keys(models)); by = String)
        m = models[site]
        n_active = m.sparsity_level
        r2 = haskey(r2_by_site, site) ? _num_to_latex(Float64(r2_by_site[site])) : "-"
        rn = haskey(residual_by_site, site) ? _num_to_latex(Float64(residual_by_site[site])) : _num_to_latex(m.residual_norm)
        push!(lines, string(_latex_symbol(site), " & ", n_active, " & ", r2, " & ", rn, " ", "\\\\"))
    end

    push!(lines, bot_rule)
    push!(lines, "\\end{tabular}")
    return join(lines, "\n")
end