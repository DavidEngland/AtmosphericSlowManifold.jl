# src/Discovery/ConstraintBuilder.jl
abstract type AbstractPhysicalConstraint end

struct PositivityConstraint <: AbstractPhysicalConstraint
    target_feature::AbstractBasisFeature
end

struct MonotonicityConstraint <: AbstractPhysicalConstraint
    target_feature::AbstractBasisFeature
    respect_to::AbstractBasisFeature
end

struct EnergyConstraint <: AbstractPhysicalConstraint
    dissipation_term::AbstractBasisFeature
end

struct PhysicalConstraintMatrix
    A_ineq::Matrix{Float64}
    b_ineq::Vector{Float64}
end

function _feature_index(library::Vector{AbstractBasisFeature}, target::AbstractBasisFeature)
    idx = findfirst(f -> typeof(f) == typeof(target) && f == target, library)
    idx === nothing && throw(ArgumentError("Feature not found in library: $(target)"))
    return idx
end

function assemble_constraint_matrix(
    constraints::Vector{<:AbstractPhysicalConstraint},
    library::Vector{AbstractBasisFeature},
    eval_grid::Matrix{Float64},
)
    n_features = length(library)
    size(eval_grid, 2) == n_features || throw(ArgumentError("eval_grid columns must match library size"))

    rows = Vector{Vector{Float64}}()
    rhs = Float64[]

    for c in constraints
        if c isa PositivityConstraint
            i = _feature_index(library, c.target_feature)
            row = zeros(n_features)
            row[i] = 1.0
            push!(rows, row)
            push!(rhs, 0.0)
        elseif c isa MonotonicityConstraint
            i = _feature_index(library, c.target_feature)
            j = _feature_index(library, c.respect_to)
            sign_proxy = sign(sum(eval_grid[:, j]))
            sign_proxy == 0.0 && (sign_proxy = 1.0)
            row = zeros(n_features)
            row[i] = sign_proxy
            push!(rows, row)
            push!(rhs, 0.0)
        elseif c isa EnergyConstraint
            i = _feature_index(library, c.dissipation_term)
            row = zeros(n_features)
            row[i] = 1.0
            push!(rows, row)
            push!(rhs, 0.0)
        else
            throw(ArgumentError("Unsupported physical constraint type: $(typeof(c))"))
        end
    end

    if isempty(rows)
        return PhysicalConstraintMatrix(zeros(0, n_features), Float64[])
    end

    A = reduce(vcat, permutedims.(rows))
    b = collect(rhs)
    return PhysicalConstraintMatrix(A, b)
end
