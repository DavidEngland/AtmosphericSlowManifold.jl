module AtmosphericSlowManifold

using ModelingToolkit
using Symbolics
using MethodOfLines
using DifferentialEquations
using DomainSets
using LinearAlgebra
using JuMP

include("Manifold/ManifoldState.jl")
include("Manifold/GSPTDiagnostics.jl")
include("Geometry/Geometry.jl")

include("Closures/Interface.jl")
include("Closures/WSINDyClosure.jl")
include("Closures/MOSTClosure.jl")

include("Observation/DataIngestion.jl")
include("Observation/SpectralBLTransform.jl")

include("System/SurfaceBoundary.jl")
include("System/PrognosticPDE.jl")

include("Discretization/Interface.jl")
include("Discretization/StretchedGrid.jl")
include("Discretization/Backends/MethodOfLinesFD.jl")
include("Discretization/Backends/SpectralBLGalerkin.jl")

include("Discovery/Discovery.jl")
include("Calibration/HierarchicalTuring.jl")

export ManifoldState, FoldConstraint, fold_residual, fold_transversality
export Geometry

export AbstractClosure, WSINDyClosure, MOSTClosure
export eddy_momentum, eddy_heat, surface_flux

export AbstractDiscretization, MethodOfLinesFD, SpectralBLGalerkin
export generate_stretched_grid, solve_scm
export SpectralNonlinearTensors, precompute_nonlinear_tensors

export ObservationTable, read_tower_csv, read_tower_netcdf, project_to_gegenbauer

export GegenbauerBasis, build_weak_library, fit_wsindy_jump, extract_closure
export discover_closure
export discover
export AbstractBasisFeature, StateVariable, SpatialDerivative, DiagnosticVariable
export BasisOperator, OperatorTerm, DiscoveredModel, to_mtk_expression, get_feature_symbolic
export FeatureLibrary, build_feature_library
export AbstractPhysicalConstraint, PositivityConstraint, MonotonicityConstraint, EnergyConstraint
export PhysicalConstraintMatrix, assemble_constraint_matrix
export AbstractTestFunctionFamily, GegenbauerFamily, BSplineFamily
export evaluate_test_function, evaluate_dt_test_function, evaluate_dz2_test_function
export WeakFormMatrix, assemble_weak_system
export AbstractSparseOptimizer, STRidge, ConstrainedQP, solve_sparse_regression

export build_pde_system, default_surface_flux

export verify_closure

end # module AtmosphericSlowManifold
