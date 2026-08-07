# src/Calibration/HierarchicalTuring.jl
using LinearAlgebra

export CalibrationConfig, HierarchicalCalibrationResult, calibrate_hierarchical

struct CalibrationConfig
    global_scale::Float64
    site_scale::Float64
    samples::Int
    chains::Int
    target_accept::Float64
end

function CalibrationConfig(;
    global_scale::Float64 = 1.0,
    site_scale::Float64 = 0.2,
    samples::Int = 1000,
    chains::Int = 4,
    target_accept::Float64 = 0.8,
)
    return CalibrationConfig(global_scale, site_scale, samples, chains, target_accept)
end

struct HierarchicalCalibrationResult
    global_mean::Vector{Float64}
    global_std::Vector{Float64}
    site_means::Dict{Symbol, Vector{Float64}}
    site_stds::Dict{Symbol, Vector{Float64}}
    diagnostics::Dict{Symbol, Any}
end

"""
Calibrate a DiscoveredModel across multiple observation sites simultaneously
using hierarchical Bayesian inference in Turing.jl.
"""
function calibrate_hierarchical(
    model::DiscoveredModel{Float64},
    site_observations::Dict{Symbol, ObservationTable},
    config::CalibrationConfig = CalibrationConfig(),
)
    Base.find_package("Turing") !== nothing ||
        throw(ErrorException("Turing.jl is required for hierarchical calibration."))

    Core.eval(@__MODULE__, :(import Turing, Distributions))

    site_names = collect(keys(site_observations))
    isempty(site_names) && throw(ArgumentError("site_observations dictionary cannot be empty."))

    X_sites = Matrix{Float64}[]
    y_sites = Vector{Float64}[]

    for site in site_names
        obs = site_observations[site]
        ysym = _cal_target_symbol(model, obs)
        push!(y_sites, Float64.(obs.columns[ysym]))
        push!(X_sites, _cal_design_matrix(obs, model))
    end

    n_terms = length(model.terms)
    n_sites = length(site_names)

    if n_terms == 0
        return HierarchicalCalibrationResult(
            Float64[],
            Float64[],
            Dict(s => Float64[] for s in site_names),
            Dict(s => Float64[] for s in site_names),
            Dict(:status => :ok, :n_sites => n_sites, :n_terms => 0),
        )
    end

    kernel = Core.eval(
        @__MODULE__,
        quote
            Turing.@model function hierarchical_pde_kernel(X_list, y_list, p, S, cfg)
                θ_global ~ Distributions.MvNormal(zeros(p), (cfg.global_scale^2) * I)
                σ_site ~ Distributions.Exponential(cfg.site_scale)
                σ_obs ~ Distributions.Gamma(2.0, 0.1)

                θ_sites = Vector{Vector{Real}}(undef, S)
                for s in 1:S
                    θ_sites[s] ~ Distributions.MvNormal(θ_global, (σ_site^2) * I)
                    y_pred = X_list[s] * θ_sites[s]
                    y_list[s] ~ Distributions.MvNormal(y_pred, (σ_obs^2) * I)
                end
            end
            hierarchical_pde_kernel
        end,
    )

    prob_model = kernel(X_sites, y_sites, n_terms, n_sites, config)
    sampler = Turing.NUTS(config.target_accept)
    chain = Turing.sample(prob_model, sampler, Turing.MCMCSerial(), config.samples, config.chains)

    draws_mat = _draw_matrix_from_chain(chain)
    names_vec = _chain_param_names(chain)

    global_mean = zeros(n_terms)
    global_std = zeros(n_terms)
    site_means = Dict{Symbol, Vector{Float64}}()
    site_stds = Dict{Symbol, Vector{Float64}}()

    for s in site_names
        site_means[s] = zeros(n_terms)
        site_stds[s] = zeros(n_terms)
    end

    for j in 1:n_terms
        g_idx = findfirst(n -> occursin("θ_global[$j]", String(n)), names_vec)
        if g_idx !== nothing
            vals = vec(draws_mat[g_idx, :])
            global_mean[j] = sum(vals) / length(vals)
            global_std[j] = sqrt(max(sum((vals .- global_mean[j]) .^ 2) / (length(vals) - 1), 0.0))
        end

        for (s_idx, s_name) in enumerate(site_names)
            s_target = "θ_sites[$s_idx][$j]"
            st_idx = findfirst(n -> occursin(s_target, String(n)), names_vec)
            if st_idx !== nothing
                vals = vec(draws_mat[st_idx, :])
                m = sum(vals) / length(vals)
                site_means[s_name][j] = m
                site_stds[s_name][j] = sqrt(max(sum((vals .- m) .^ 2) / (length(vals) - 1), 0.0))
            end
        end
    end

    diagnostics = Dict{Symbol, Any}(
        :status => :ok,
        :backend => :HierarchicalTuring,
        :n_sites => n_sites,
        :n_terms => n_terms,
        :chain => chain,
        :config => config,
    )

    return HierarchicalCalibrationResult(global_mean, global_std, site_means, site_stds, diagnostics)
end
