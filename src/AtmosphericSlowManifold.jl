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

include("Discovery/SymbolicVerification.jl")
include("Calibration/HierarchicalTuring.jl")

export ManifoldState, FoldConstraint, fold_residual, fold_transversality

export AbstractClosure, WSINDyClosure, MOSTClosure
export eddy_momentum, eddy_heat, surface_flux

export AbstractDiscretization, MethodOfLinesFD, SpectralBLGalerkin
export generate_stretched_grid, solve_scm

export ObservationTable, read_tower_csv, read_tower_netcdf, project_to_gegenbauer

export build_pde_system, default_surface_flux

export verify_closure

end # module AtmosphericSlowManifold
