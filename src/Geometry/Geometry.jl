module Geometry

using LinearAlgebra

include("Jacobians.jl")
using .Jacobians: JacobianModel, JacobianCache, evaluate_jacobian, evaluate_det, evaluate_adjugate, compute_fast_jacobian, compute_adjugate, evaluate_tangent_space

include("CriticalManifold.jl")
include("DesingularizedFlow.jl")
include("Fenichel.jl")
include("FoldTracking.jl")
include("CanardDetection.jl")
include("Continuation.jl")

using .Continuation: AbstractContinuationAlgorithm, PseudoArclength, ContinuationBranch, continue_manifold, continue_set

export ManifoldPoint, CriticalManifoldSurface, find_manifold_point, solve_critical_surface
export AbstractInvariantSet
export fold_indicator, track_fold_curve
export desingularized_vector_field, slow_flow_vector
export find_desingularized_singular_points
export FoldedSingularity, SingularType, CanardSegment, classify_folded_singularity, classify_singular_type, detect_folded_singularity, build_canard_segment
export HyperbolicityReport, fenichel_metrics, hyperbolicity_profile
export JacobianModel, JacobianCache, evaluate_jacobian, evaluate_det, evaluate_adjugate
export compute_fast_jacobian, compute_adjugate, evaluate_tangent_space
export AbstractContinuationAlgorithm, PseudoArclength, ContinuationBranch, continue_manifold, continue_set
export FoldCurve, track_fold_line
export Jacobians, Continuation

function finite_difference_jacobian_y(f_fast, x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; eps::Float64 = 1e-6, params = Dict{Symbol, Float64}())
    return compute_fast_jacobian(f_fast, x, y; eps = eps, params = params)
end

end # module Geometry
