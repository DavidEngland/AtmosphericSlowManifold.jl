module Calibration

using ..AtmosphericSlowManifold: DiscoveredModel, ObservationTable
using ..AtmosphericSlowManifold: BasisOperator, OperatorTerm, AbstractBasisFeature
using ..AtmosphericSlowManifold: StateVariable, DiagnosticVariable, SpatialDerivative
using LinearAlgebra

include("Interface.jl")
include("Backends/MaximumLikelihood.jl")
include("Backends/BayesianMCMC.jl")
include("Backends/VariationalInference.jl")

end # module Calibration
