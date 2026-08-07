module CampaignDiagnosticsModule

export CampaignDiagnostics

struct CampaignDiagnostics{T<:AbstractFloat}
    campaign_name::String
    num_samples::Int
    rmse_u::T
    rmse_theta::T
    r2_closure::T
    mean_closure_residual::T
    energy_imbalance_rms::T
    integrated_dissipation::T
    transversality_mean::T
    transversality_std::T
    manifold_distance_mean::T
    fold_distance_min::T
end

end # module CampaignDiagnosticsModule
