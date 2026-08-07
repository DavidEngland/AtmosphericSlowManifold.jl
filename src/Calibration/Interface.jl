# src/Calibration/Interface.jl
export AbstractCalibrationAlgorithm
export BayesianMCMC, MaximumLikelihood, VariationalInference
export CalibrationResult, calibrate, dispatch_calibrate

abstract type AbstractCalibrationAlgorithm end

struct BayesianMCMC <: AbstractCalibrationAlgorithm
    samples::Int
    chains::Int
    target_accept::Float64
    sampler_type::Symbol
end

struct MaximumLikelihood <: AbstractCalibrationAlgorithm
    optimizer::Symbol
    max_iter::Int
    tol::Float64
end

struct VariationalInference <: AbstractCalibrationAlgorithm
    samples::Int
    max_iter::Int
    optimizer::Symbol
end

struct CalibrationResult{A<:AbstractCalibrationAlgorithm}
    model::DiscoveredModel{Float64}
    parameters::Dict{Symbol, Vector{Float64}}
    diagnostics::Dict{Symbol, Any}
    algorithm::A
end

function calibrate(
    model::DiscoveredModel,
    obs::ObservationTable;
    algorithm::AbstractCalibrationAlgorithm = BayesianMCMC(1000, 4, 0.8, :NUTS),
)
    return dispatch_calibrate(model, obs, algorithm)
end

function dispatch_calibrate(
    model::DiscoveredModel,
    obs::ObservationTable,
    alg::AbstractCalibrationAlgorithm,
)
    error("Calibration backend $(typeof(alg)) is not implemented.")
end
