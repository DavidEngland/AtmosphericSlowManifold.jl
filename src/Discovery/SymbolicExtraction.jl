abstract type AbstractBasisFeature end

struct StateVariable <: AbstractBasisFeature
    name::Symbol
end

struct SpatialDerivative <: AbstractBasisFeature
    variable::Symbol
    order::Int
end

struct DiagnosticVariable <: AbstractBasisFeature
    name::Symbol
end

struct BasisOperator
    feature::AbstractBasisFeature
    power::Float64
end

struct OperatorTerm{T}
    coefficient::T
    basis::Vector{BasisOperator}
end

struct DiscoveredModel{T}
    target_variable::Symbol
    terms::Vector{OperatorTerm{T}}
    residual_norm::Float64
    sparsity_level::Int
end

function _coefficient_to_num(c)
    if c isa Num
        return c
    elseif c isa Number
        return Num(c)
    else
        return Num(float(c))
    end
end

function get_feature_symbolic(f::StateVariable, var_map::Dict{Symbol, Num})
    haskey(var_map, f.name) || throw(KeyError(f.name))
    return var_map[f.name]
end

function get_feature_symbolic(f::DiagnosticVariable, var_map::Dict{Symbol, Num})
    haskey(var_map, f.name) || throw(KeyError(f.name))
    return var_map[f.name]
end

function get_feature_symbolic(f::SpatialDerivative, var_map::Dict{Symbol, Num})
    key = Symbol("d$(f.order)_$(f.variable)")
    haskey(var_map, key) || throw(KeyError(key))
    return var_map[key]
end

function to_mtk_expression(term::OperatorTerm{T}, var_map::Dict{Symbol, Num}) where {T}
    expr = _coefficient_to_num(term.coefficient)
    for b in term.basis
        base_num = get_feature_symbolic(b.feature, var_map)
        expr *= (b.power == 1.0) ? base_num : base_num^b.power
    end
    return expr
end

function to_mtk_expression(model::DiscoveredModel{T}, var_map::Dict{Symbol, Num}) where {T}
    expr = Num(0.0)
    for term in model.terms
        expr += to_mtk_expression(term, var_map)
    end
    return expr
end
