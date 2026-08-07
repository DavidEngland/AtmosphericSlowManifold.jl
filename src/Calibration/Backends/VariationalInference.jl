# src/Calibration/Backends/VariationalInference.jl
function dispatch_calibrate(
    model::DiscoveredModel{Float64},
    obs::ObservationTable,
    alg::VariationalInference,
)
    ysym = _cal_target_symbol(model, obs)
    y = Float64.(obs.columns[ysym])
    X = _cal_design_matrix(obs, model)
    p = size(X, 2)

    if p == 0
        params = Dict{Symbol, Vector{Float64}}(
            :variational_mean => Float64[],
            :variational_scale => Float64[],
        )
        diagnostics = Dict{Symbol, Any}(
            :status => :ok,
            :backend => :VariationalInference,
            :optimizer => alg.optimizer,
            :max_iter => alg.max_iter,
            :target_symbol => ysym,
            :n_obs => length(y),
            :n_terms => 0,
            :converged => true,
            :final_elbo => 0.0,
            :elbo_trace => Float64[],
        )
        return CalibrationResult{VariationalInference}(model, params, diagnostics, alg)
    end

    # Mean-field Gaussian q(θ) = N(μ, diag(σ²)) with fixed Gaussian prior θ ~ N(0, τ²I).
    noise_var = max(sum((y .- (sum(y) / max(length(y), 1))) .^ 2) / max(length(y), 1), 1e-8)
    tau2 = 1.0
    gram_diag = vec(sum(X .^ 2; dims = 1))

    μ = _cal_fit_ridge(X, y, max(1e-10, sqrt(noise_var) * 1e-6))
    logσ = fill(log(0.1), p)

    elbo_trace = Float64[]
    converged = false

    # Optimizer state (ADAM default).
    m_μ = zeros(p)
    v_μ = zeros(p)
    m_s = zeros(p)
    v_s = zeros(p)
    β1 = 0.9
    β2 = 0.999
    eps = 1e-8
    lr = alg.optimizer == :SGD ? 5e-2 : 1e-2

    for it in 1:max(1, alg.max_iter)
        σ2 = exp.(2.0 .* logσ)

        resid = y .- (X * μ)
        res2 = dot(resid, resid)

        # ELBO up to additive constants.
        expected_sq_err = res2 + dot(gram_diag, σ2)
        kl = 0.5 * sum((σ2 .+ μ .^ 2) ./ tau2 .- 1.0 .- log.(σ2 ./ tau2))
        elbo = -0.5 * expected_sq_err / noise_var - kl
        push!(elbo_trace, elbo)

        # Gradients for minimizing negative ELBO.
        gμ = -(X' * resid) / noise_var + μ / tau2
        gs = (gram_diag .* σ2) / noise_var + (σ2 / tau2) .- 1.0

        gnorm = sqrt(sum(gμ .^ 2) + sum(gs .^ 2))
        if gnorm < 1e-6
            converged = true
            break
        end

        if alg.optimizer == :SGD
            μ .-= lr .* gμ
            logσ .-= lr .* gs
        else
            m_μ .= β1 .* m_μ .+ (1.0 - β1) .* gμ
            v_μ .= β2 .* v_μ .+ (1.0 - β2) .* (gμ .^ 2)
            m_s .= β1 .* m_s .+ (1.0 - β1) .* gs
            v_s .= β2 .* v_s .+ (1.0 - β2) .* (gs .^ 2)

            mhat_μ = m_μ ./ (1.0 - β1^it)
            vhat_μ = v_μ ./ (1.0 - β2^it)
            mhat_s = m_s ./ (1.0 - β1^it)
            vhat_s = v_s ./ (1.0 - β2^it)

            μ .-= lr .* mhat_μ ./ (sqrt.(vhat_μ) .+ eps)
            logσ .-= lr .* mhat_s ./ (sqrt.(vhat_s) .+ eps)
        end

        # Keep scales numerically stable.
        logσ .= clamp.(logσ, log(1e-6), log(10.0))
    end

    σ = exp.(logσ)
    posterior_mean = copy(μ)
    posterior_scale = copy(σ)

    yhat = X * posterior_mean
    residual = yhat .- y
    calibrated_model = _cal_update_model(model, posterior_mean, norm(residual))

    params = Dict{Symbol, Vector{Float64}}(
        :variational_mean => posterior_mean,
        :variational_scale => posterior_scale,
    )
    ndraws = max(1, alg.samples)
    for i in 1:p
        params[Symbol("c_", i)] = posterior_mean[i] .+ posterior_scale[i] .* randn(ndraws)
    end

    diagnostics = Dict{Symbol, Any}(
        :status => :ok,
        :backend => :VariationalInference,
        :optimizer => alg.optimizer,
        :max_iter => alg.max_iter,
        :target_symbol => ysym,
        :n_obs => length(y),
        :n_terms => p,
        :converged => converged,
        :final_elbo => isempty(elbo_trace) ? 0.0 : elbo_trace[end],
        :elbo_trace => elbo_trace,
        :noise_var => noise_var,
    )

    return CalibrationResult{VariationalInference}(calibrated_model, params, diagnostics, alg)
end
