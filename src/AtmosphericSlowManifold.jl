module AtmosphericSlowManifold

using ModelingToolkit
using Symbolics
using MethodOfLines
using DifferentialEquations
using DomainSets
using LinearAlgebra
using JuMP
using JSON3
using CSV
using DataFrames
using Dates

include("Manifold/ManifoldState.jl")
include("Manifold/GSPTDiagnostics.jl")
include("Geometry/Geometry.jl")

include("Closures/Interface.jl")
include("Closures/SmoothOperators.jl")
include("Closures/WSINDyClosure.jl")
include("Closures/MOSTClosure.jl")
include("Closures/PhysicalSimilarityClosure.jl")
include("Closures/Z0HRClosure.jl")

include("Observation/DataIngestion.jl")
include("Observation/SpectralBLTransform.jl")

include("Forcing/SurfaceForcing.jl")
include("System/SurfaceBoundary.jl")
include("System/PrognosticPDE.jl")
include("Diagnostics/Diagnostics.jl")

include("Discretization/Interface.jl")
include("Discretization/StretchedGrid.jl")
include("Discretization/Backends/MethodOfLinesFD.jl")
include("Discretization/Backends/SpectralBLGalerkin.jl")

include("Discovery/Discovery.jl")
include("Calibration/Calibration.jl")
include("System/ExportUtilities.jl")

using .Calibration: AbstractCalibrationAlgorithm, BayesianMCMC, MaximumLikelihood, VariationalInference
using .Calibration: CalibrationResult, calibrate, dispatch_calibrate
using .Calibration: CalibrationConfig, HierarchicalCalibrationResult, calibrate_hierarchical
using .Calibration: evaluate_profile_uncertainty
using .Geometry: FoldedSingularity, classify_folded_singularity, detect_canard_trajectories
using .ExportUtilities: export_to_csv, export_to_json, export_to_netcdf

export ManifoldState, FoldConstraint, fold_residual, fold_transversality
export Geometry
export FoldedSingularity, classify_folded_singularity, detect_canard_trajectories

export AbstractClosure, AbstractAtmosphericClosure, WSINDyClosure, MOSTClosure, PhysicalSimilarityClosure, Z0HRClosure
export eddy_momentum, eddy_heat, surface_flux
export evaluate_diffusivity_profile!, evaluate_heat_diffusivity_profile!
export smooth_max, smooth_min, smooth_floor, z0hr_stability_functions

export AbstractDiscretization, MethodOfLinesFD, SpectralBLGalerkin
export generate_stretched_grid, solve_scm
export SpectralNonlinearTensors, precompute_nonlinear_tensors
export ModalBudgetDiagnostic, evaluate_modal_budget
export SurfaceForcing, load_surface_forcing, interp_forcing

export ObservationTable, read_tower_csv, read_tower_netcdf, read_observation_data, project_to_gegenbauer
export resolve_sibling_data_dir, find_data_files

export GegenbauerBasis, build_weak_library, fit_wsindy_jump, extract_closure
export constrained_stlsq
export discover_closure
export discover
export AbstractBasisFeature, StateVariable, SpatialDerivative, DiagnosticVariable
export BasisOperator, OperatorTerm, DiscoveredModel, to_mtk_expression, get_feature_symbolic
export to_latex, latex_term_table, latex_site_summary_table
export aic, bic, model_aic, model_bic, compute_pareto_front, kfold_cv_residual
export FeatureLibrary, build_feature_library
export SyntheticBenchmarkData
export manufactured_linear_solution, manufactured_nonlinear_solution
export linear_diffusivity, nonlinear_diffusivity
export generate_manufactured_field
export additive_gaussian_noise, ar1_temporal_noise, multiplicative_sensor_noise, apply_missing_data
export coefficient_l2_error, precision_recall, structural_hamming_distance, false_discovery_rate, equation_sparsity, snr_db, evaluate_recovery_metrics
export IdentifiabilityReport, compute_gram_matrix, compute_mutual_coherence, compute_condition_number, compute_vif, compute_parameter_covariance, analyze_identifiability, prune_by_mutual_coherence
export CampaignData, LOSOSplit, LOSOResult, LOSOSummary
export create_loso_splits, evaluate_out_of_sample, run_loso_cross_validation
export load_campaign_data, load_all_campaigns, run_artifact_loso, export_loso_table
export AbstractPhysicalConstraint, PositivityConstraint, MonotonicityConstraint, EnergyConstraint
export PhysicalConstraintMatrix, assemble_constraint_matrix
export AbstractTestFunctionFamily, GegenbauerFamily, BSplineFamily
export evaluate_test_function, evaluate_dt_test_function, evaluate_dz2_test_function
export WeakFormMatrix, assemble_weak_system
export AbstractSparseOptimizer, STRidge, ConstrainedQP, solve_sparse_regression

export Calibration
export AbstractCalibrationAlgorithm, BayesianMCMC, MaximumLikelihood, VariationalInference
export CalibrationResult, calibrate, dispatch_calibrate
export CalibrationConfig, HierarchicalCalibrationResult, calibrate_hierarchical
export evaluate_profile_uncertainty

export build_pde_system, default_surface_flux, default_surface_heat_flux, surface_boundary_conditions
export export_to_csv, export_to_json, export_to_netcdf
export Diagnostics

export verify_closure
export Z0HRClosure, z0hr_stability_functions, eddy_momentum, eddy_heat, surface_flux, evaluate_diffusivity_profile!, evaluate_heat_diffusivity_profile!

end # module AtmosphericSlowManifold