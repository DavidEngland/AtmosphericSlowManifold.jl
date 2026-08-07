# src/Discovery/LibraryBuilder.jl
struct FeatureLibrary
    features::Vector{AbstractBasisFeature}
    evaluators::Vector{Function}
end

function _feature_key(f::StateVariable)
    return f.name
end

function _feature_key(f::DiagnosticVariable)
    return f.name
end

function _feature_key(f::SpatialDerivative)
    return Symbol("d$(f.order)_$(f.variable)")
end

function build_feature_library(
    states::Vector{Symbol},
    diagnostics::Vector{Symbol},
    max_derivative_order::Int,
)
    max_derivative_order >= 0 || throw(ArgumentError("max_derivative_order must be non-negative"))

    features = AbstractBasisFeature[]

    for s in states
        push!(features, StateVariable(s))
        for order in 1:max_derivative_order
            push!(features, SpatialDerivative(s, order))
        end
    end

    for d in diagnostics
        push!(features, DiagnosticVariable(d))
    end

    evaluators = Function[]
    for f in features
        key = _feature_key(f)
        push!(evaluators, sample -> begin
            haskey(sample, key) || throw(KeyError(key))
            return float(sample[key])
        end)
    end

    return FeatureLibrary(features, evaluators)
end
