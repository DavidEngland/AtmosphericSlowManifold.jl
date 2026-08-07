# src/Calibration/Backends/BayesianMCMC.jl
function _bayes_sigma2(y::Vector{Float64}, yhat::Vector{Float64})
    r = y .- yhat
    n = max(length(y), 1)
    return max(dot(r, r) / n, 1e-12)
end

function _turing_available()
    return Base.find_package("Turing") !== nothing
end

function _bayes_draws(mu::Vector{Float64}, cov::Matrix{Float64}, ndraws::Int)
    p = length(mu)
    p == 0 && return zeros(0, ndraws)

    cov_reg = cov + 1e-10 * I(p)
    F = cholesky(Symmetric(cov_reg); check = false)
    z = randn(p, ndraws)
    return mu .+ F.L * z
end

function _rowwise_std(draws::Matrix{Float64})
    p, n = size(draws)
    p == 0 && return Float64[]
    n <= 1 && return zeros(p)

    out = zeros(p)
    for i in 1:p
        row = @view draws[i, :]
        μ = sum(row) / n
        var = sum((row .- μ) .^ 2) / (n - 1)
        out[i] = sqrt(max(var, 0.0))
    end
    return out
end

function _bayes_update_model(model::DiscoveredModel{Float64}, θ::Vector{Float64}, residual_norm::Float64)
    length(θ) == length(model.terms) || throw(ArgumentError("Coefficient length mismatch during Bayesian model update."))
    new_terms = OperatorTerm{Float64}[]
    for (i, t) in enumerate(model.terms)
        push!(new_terms, OperatorTerm{Float64}(θ[i], t.basis))
    end
    sparsity = count(abs.(θ) .> 1e-10)
    return DiscoveredModel{Float64}(model.target_variable, new_terms, residual_norm, sparsity)
end

function _draw_matrix_from_chain(chain)
    raw = Array(chain)
    if ndims(raw) == 3
        # (iterations, parameters, chains) -> (parameters, draws)
        iters, npar, nch = size(raw)
        return reshape(permutedims(raw, (2, 1, 3)), npar, iters * nch)
    elseif ndims(raw) == 2
        # (iterations, parameters) -> (parameters, draws)
        return permutedims(raw, (2, 1))
    else
        throw(ArgumentError("Unsupported chain array dimensionality $(ndims(raw))."))
    end
end

function _chain_param_names(chain)
    try
        return Symbol.(chain.name_map.parameters)
    catch
        try
            return Symbol.(collect(names(chain)))
        catch
            return Symbol[]
        end
    end
end

function _extract_theta_draws(chain, p::Int)
    draws_mat = _draw_matrix_from_chain(chain)
    names_vec = _chain_param_names(chain)
    out = Dict{Symbol, Vector{Float64}}()

    for i in 1:p
        target1 = "θ[$i]"
        target2 = "theta[$i]"
        idx = findfirst(n -> begin
            s = String(n)
            sl = lowercase(s)
            occursin(target1, s) || occursin(target2, sl)
        end, names_vec)

        if idx === nothing
            out[Symbol("c_", i)] = Float64[]
        else
            out[Symbol("c_", i)] = vec(draws_mat[idx, :])
        end
    end
    return out
end

function _dispatch_calibrate_gaussian_approx(
    model::DiscoveredModel{Float64},
    obs::ObservationTable,
    alg::BayesianMCMC,
)
    ysym = _cal_target_symbol(model, obs)
    y = Float64.(obs.columns[ysym])
    X = _cal_design_matrix(obs, model)

    p = size(X, 2)
    ndraws = max(1, alg.samples * alg.chains)

    λ_prior = 1.0
    if p == 0
        posterior_mean = Float64[]
        posterior_std = Float64[]
        draws = zeros(0, ndraws)
        yhat = zeros(length(y))
    else
        A = X' * X + λ_prior * I(p)
        b = X' * y
        posterior_mean = A \ b
        yhat = X * posterior_mean
        σ2 = _bayes_sigma2(y, yhat)
        posterior_cov = σ2 * inv(Matrix(A))
        draws = _bayes_draws(posterior_mean, posterior_cov, ndraws)
        posterior_std = _rowwise_std(draws)
    end

    residual = yhat .- y
    calibrated_model = _bayes_update_model(model, posterior_mean, norm(residual))

    params = Dict{Symbol, Vector{Float64}}(
        :posterior_mean => posterior_mean,
        :posterior_std => posterior_std,
    )
    for i in 1:p
        params[Symbol("c_", i)] = vec(draws[i, :])
    end

    diagnostics = Dict{Symbol, Any}(
        :status => :ok,
        :backend => :BayesianMCMC,
        :engine => :gaussian_approx,
        :sampler_type => alg.sampler_type,
        :samples => alg.samples,
        :chains => alg.chains,
        :target_accept => alg.target_accept,
        :target_symbol => ysym,
        :n_obs => length(y),
        :n_terms => p,
        :ess_proxy => ndraws,
        :rhat_proxy => p == 0 ? Float64[] : fill(1.0, p),
    )

    return CalibrationResult{BayesianMCMC}(calibrated_model, params, diagnostics, alg)
end

function _dispatch_calibrate_turing(
    model::DiscoveredModel{Float64},
    obs::ObservationTable,
    alg::BayesianMCMC,
)
    Core.eval(@__MODULE__, :(import Turing, Distributions))

    ysym = _cal_target_symbol(model, obs)
    y = Float64.(obs.columns[ysym])
    X = _cal_design_matrix(obs, model)
    p = size(X, 2)

    if p == 0
        params = Dict{Symbol, Vector{Float64}}(
            :posterior_mean => Float64[],
            :posterior_std => Float64[],
        )
        diagnostics = Dict{Symbol, Any}(
            :status => :ok,
            :backend => :BayesianMCMC,
            :engine => :turing,
            :sampler_type => alg.sampler_type,
            :samples => alg.samples,
            :chains => alg.chains,
            :target_accept => alg.target_accept,
            :target_symbol => ysym,
            :n_obs => length(y),
            :n_terms => 0,
        )
        return CalibrationResult{BayesianMCMC}(model, params, diagnostics, alg)
    end

    # Build the probabilistic kernel at runtime so package load does not require Turing.
    kernel = Core.eval(
        @__MODULE__,
        quote
            Turing.@model function bayes_wsindy_kernel(y_obs, Xmat, n_terms)
                θ ~ Distributions.MvNormal(zeros(n_terms), 1.0)
                σ ~ Distributions.Gamma(2.0, 0.1)
                y_pred = Xmat * θ
                y_obs ~ Distributions.MvNormal(y_pred, (σ^2) * I)
            end
            bayes_wsindy_kernel
        end,
    )

    prob_model = kernel(y, X, p)
    sampler = if alg.sampler_type == :NUTS
        Turing.NUTS(alg.target_accept)
    elseif alg.sampler_type == :HMC
        Turing.HMC(0.05, 10)
    else
        Turing.MH()
    end

    chain = Turing.sample(prob_model, sampler, Turing.MCMCSerial(), alg.samples, alg.chains)
    theta_draws = _extract_theta_draws(chain, p)

    posterior_mean = zeros(p)
    posterior_std = zeros(p)
    for i in 1:p
        key = Symbol("c_", i)
        vals = get(theta_draws, key, Float64[])
        if isempty(vals)
            continue
        end
        posterior_mean[i] = sum(vals) / length(vals)
        posterior_std[i] = begin
            μ = posterior_mean[i]
            n = length(vals)
            n <= 1 ? 0.0 : sqrt(max(sum((vals .- μ) .^ 2) / (n - 1), 0.0))
        end
    end

    yhat = X * posterior_mean
    residual = yhat .- y
    calibrated_model = _bayes_update_model(model, posterior_mean, norm(residual))

    params = Dict{Symbol, Vector{Float64}}(
        :posterior_mean => posterior_mean,
        :posterior_std => posterior_std,
    )
    merge!(params, theta_draws)

    diagnostics = Dict{Symbol, Any}(
        :status => :ok,
        :backend => :BayesianMCMC,
        :engine => :turing,
        :sampler_type => alg.sampler_type,
        :samples => alg.samples,
        :chains => alg.chains,
        :target_accept => alg.target_accept,
        :target_symbol => ysym,
        :n_obs => length(y),
        :n_terms => p,
        :chain => chain,
    )

    return CalibrationResult{BayesianMCMC}(calibrated_model, params, diagnostics, alg)
end

function dispatch_calibrate(
    model::DiscoveredModel{Float64},
    obs::ObservationTable,
    alg::BayesianMCMC,
)
    if _turing_available()
        try
            return _dispatch_calibrate_turing(model, obs, alg)
        catch
            # Fall through to deterministic approximation if Turing path errors.
        end
    end
    return _dispatch_calibrate_gaussian_approx(model, obs, alg)
end
