struct FoldConstraint
    expr::Num
end

"""Evaluate fold expression with substitutions and return a numeric residual."""
function fold_residual(fold::FoldConstraint, substitutions::Dict{Num, <:Real})
    return Symbolics.value(Symbolics.substitute(fold.expr, substitutions))
end

"""Differentiate the fold expression with respect to a symbolic variable."""
function fold_transversality(fold::FoldConstraint, variable::Num)
    return Symbolics.derivative(fold.expr, variable)
end
