# src/Discovery/Discovery.jl
include("SymbolicExtraction.jl")
include("LibraryBuilder.jl")
include("ConstraintBuilder.jl")
include("TestFunctions.jl")
include("WeakForms.jl")
include("SparseRegression.jl")
include("WSINDyEngine.jl")
include("SyntheticBenchmarks.jl")
include("Identifiability.jl")
include("LOSOValidation.jl")
include("SymbolicVerification.jl")
include("LaTeXExporter.jl")
include("ModelSelection.jl")

export CampaignData, LOSOSplit, LOSOResult, LOSOSummary
export create_loso_splits, evaluate_out_of_sample, run_loso_cross_validation
export load_campaign_data, load_all_campaigns, run_artifact_loso
