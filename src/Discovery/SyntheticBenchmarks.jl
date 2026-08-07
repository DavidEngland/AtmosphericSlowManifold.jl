using Random
using Statistics

struct SyntheticBenchmarkData{T<:Real}
    z::Vector{T}
    t::Vector{T}
    field::Matrix{T}
    model::Symbol
end

function manufactured_linear_solution(z::Real, t::Real; lambda::Float64 = 0.25)
    return exp(-lambda * t) * sin(pi * z)
end

function manufactured_nonlinear_solution(
    z::Real,
    t::Real;
    U0::Float64 = 1.0,
    z0::Float64 = 0.05,
    p::Float64 = 0.2,
    gamma::Float64 = 0.1,
)
    z_eff = max(float(z), z0)
    return U0 * (z_eff / z0)^p * exp(-gamma * t)
end

linear_diffusivity(z::Real; nu0::Float64 = 0.05, nu1::Float64 = 0.02) = nu0 + nu1 * float(z)

function nonlinear_diffusivity(
    z::Real,
    uz::Real;
    kappa::Float64 = 0.4,
    alpha::Float64 = 5.0,
    ri::Float64 = 0.1,
)
    return (kappa^2 * float(z)^2 * float(uz)) / (1.0 + alpha * ri)
end

function generate_manufactured_field(
    z::AbstractVector{<:Real},
    t::AbstractVector{<:Real};
    model::Symbol = :linear,
    kwargs...
)
    zf = Float64.(z)
    tf = Float64.(t)
    field = zeros(Float64, length(zf), length(tf))

    if model == :linear
        for i in eachindex(zf), j in eachindex(tf)
            field[i, j] = manufactured_linear_solution(zf[i], tf[j]; kwargs...)
        end
    elseif model == :nonlinear
        for i in eachindex(zf), j in eachindex(tf)
            field[i, j] = manufactured_nonlinear_solution(zf[i], tf[j]; kwargs...)
        end
    else
        throw(ArgumentError("Unsupported manufactured model: $(model)."))
    end

    return SyntheticBenchmarkData(zf, tf, field, model)
end

function additive_gaussian_noise(
    field::AbstractMatrix{<:Real};
    sigma_rel::Float64 = 0.05,
    rng::AbstractRNG = Random.default_rng(),
)
    sigma_rel >= 0 || throw(ArgumentError("sigma_rel must be non-negative."))
    base = Float64.(field)
    sigma = sigma_rel * max(std(vec(base)), eps(Float64))
    return base .+ sigma .* randn(rng, size(base)...) 
end

function ar1_temporal_noise(
    field::AbstractMatrix{<:Real};
    sigma_rel::Float64 = 0.05,
    rho::Float64 = 0.85,
    rng::AbstractRNG = Random.default_rng(),
)
    abs(rho) < 1 || throw(ArgumentError("rho must satisfy |rho| < 1."))
    sigma_rel >= 0 || throw(ArgumentError("sigma_rel must be non-negative."))

    base = Float64.(field)
    sigma = sigma_rel * max(std(vec(base)), eps(Float64))
    noise = zeros(Float64, size(base))
    innovation_scale = sigma * sqrt(1.0 - rho^2)

    for i in axes(base, 1)
        noise[i, 1] = sigma * randn(rng)
        for j in 2:size(base, 2)
            noise[i, j] = rho * noise[i, j - 1] + innovation_scale * randn(rng)
        end
    end

    return base .+ noise
end

function multiplicative_sensor_noise(
    field::AbstractMatrix{<:Real};
    sigma_rel::Float64 = 0.05,
    rng::AbstractRNG = Random.default_rng(),
)
    sigma_rel >= 0 || throw(ArgumentError("sigma_rel must be non-negative."))
    base = Float64.(field)
    return base .* (1.0 .+ sigma_rel .* randn(rng, size(base)...))
end

function apply_missing_data(
    field::AbstractMatrix{<:Real};
    p_drop::Float64 = 0.1,
    rng::AbstractRNG = Random.default_rng(),
)
    0.0 <= p_drop <= 1.0 || throw(ArgumentError("p_drop must be in [0,1]."))
    out = Matrix{Union{Missing, Float64}}(undef, size(field)...)
    for i in eachindex(field)
        out[i] = rand(rng) < p_drop ? missing : float(field[i])
    end
    return out
end

function coefficient_l2_error(est::AbstractVector{<:Real}, truth::AbstractVector{<:Real}; eps_floor::Float64 = 1e-12)
    length(est) == length(truth) || throw(DimensionMismatch("est and truth must match in length."))
    denom = max(norm(Float64.(truth)), eps_floor)
    return norm(Float64.(est) .- Float64.(truth)) / denom
end

function precision_recall(est::AbstractVector{<:Real}, truth::AbstractVector{<:Real}; tol::Float64 = 1e-8)
    length(est) == length(truth) || throw(DimensionMismatch("est and truth must match in length."))

    est_active = abs.(Float64.(est)) .> tol
    truth_active = abs.(Float64.(truth)) .> tol

    tp = count(est_active .& truth_active)
    fp = count(est_active .& .!truth_active)
    fn = count(.!est_active .& truth_active)

    precision = tp == 0 && fp == 0 ? 1.0 : tp / max(tp + fp, 1)
    recall = tp == 0 && fn == 0 ? 1.0 : tp / max(tp + fn, 1)
    return (precision = precision, recall = recall, tp = tp, fp = fp, fn = fn)
end

function structural_hamming_distance(est::AbstractVector{<:Real}, truth::AbstractVector{<:Real}; tol::Float64 = 1e-8)
    length(est) == length(truth) || throw(DimensionMismatch("est and truth must match in length."))
    est_active = abs.(Float64.(est)) .> tol
    truth_active = abs.(Float64.(truth)) .> tol
    return count(xor.(est_active, truth_active))
end

function false_discovery_rate(est::AbstractVector{<:Real}, truth::AbstractVector{<:Real}; tol::Float64 = 1e-8)
    stats = precision_recall(est, truth; tol = tol)
    denom = stats.tp + stats.fp
    return denom == 0 ? 0.0 : stats.fp / denom
end

function equation_sparsity(coeffs::AbstractVector{<:Real}; tol::Float64 = 1e-8)
    return count(abs.(Float64.(coeffs)) .> tol)
end

function snr_db(clean::AbstractMatrix{<:Real}, noisy::AbstractMatrix{<:Real}; eps_floor::Float64 = 1e-12)
    size(clean) == size(noisy) || throw(DimensionMismatch("clean and noisy matrices must have same shape."))
    s = sum(abs2, Float64.(clean))
    n = sum(abs2, Float64.(noisy) .- Float64.(clean))
    return 10.0 * log10((s + eps_floor) / (n + eps_floor))
end

function evaluate_recovery_metrics(est::AbstractVector{<:Real}, truth::AbstractVector{<:Real}; tol::Float64 = 1e-8)
    pr = precision_recall(est, truth; tol = tol)
    return (
        l2_error = coefficient_l2_error(est, truth),
        precision = pr.precision,
        recall = pr.recall,
        shd = structural_hamming_distance(est, truth; tol = tol),
        fdr = false_discovery_rate(est, truth; tol = tol),
        sparsity = equation_sparsity(est; tol = tol),
    )
end
