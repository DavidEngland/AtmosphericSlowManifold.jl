using Statistics

export aic, bic, model_aic, model_bic, compute_pareto_front, kfold_cv_residual

"""
    aic(k, rss, n)

Compute Akaike Information Criterion:
AIC = 2k + n * log(rss / n)
"""
function aic(k::Integer, rss::Real, n::Integer)
    k >= 0 || throw(ArgumentError("k must be nonnegative."))
    n > 0 || throw(ArgumentError("n must be positive."))
    rssf = Float64(rss)
    rssf > 0 || throw(ArgumentError("rss must be positive for AIC."))
    return 2.0 * k + n * log(rssf / n)
end

"""
    bic(k, rss, n)

Compute Bayesian Information Criterion:
BIC = k * log(n) + n * log(rss / n)
"""
function bic(k::Integer, rss::Real, n::Integer)
    k >= 0 || throw(ArgumentError("k must be nonnegative."))
    n > 0 || throw(ArgumentError("n must be positive."))
    rssf = Float64(rss)
    rssf > 0 || throw(ArgumentError("rss must be positive for BIC."))
    return k * log(n) + n * log(rssf / n)
end

"""
    model_aic(model, n_obs)

Compute AIC for a discovered model using k = model sparsity and
rss = residual_norm^2 by convention.
"""
function model_aic(model::DiscoveredModel, n_obs::Integer)
    return aic(model.sparsity_level, model.residual_norm^2, n_obs)
end

"""
    model_bic(model, n_obs)

Compute BIC for a discovered model using k = model sparsity and
rss = residual_norm^2 by convention.
"""
function model_bic(model::DiscoveredModel, n_obs::Integer)
    return bic(model.sparsity_level, model.residual_norm^2, n_obs)
end

function _model_complexity(model::DiscoveredModel)
    return max(model.sparsity_level, length(model.terms))
end

function _model_objective_scores(
    models::Vector{<:DiscoveredModel};
    objective::Symbol,
    r2_values::Union{Nothing, AbstractVector{<:Real}},
)
    if objective == :rss
        return [model.residual_norm^2 for model in models], [model.residual_norm^2 for model in models]
    elseif objective == :r2
        r2_values === nothing && throw(ArgumentError("r2_values must be provided when objective=:r2."))
        length(r2_values) == length(models) || throw(ArgumentError("r2_values length must match models length."))
        scores = [-Float64(r2) for r2 in r2_values]
        raw = [Float64(r2) for r2 in r2_values]
        return scores, raw
    else
        throw(ArgumentError("Unsupported objective $(objective). Use :rss or :r2."))
    end
end

"""
    compute_pareto_front(models; objective=:rss, r2_values=nothing)

Compute a complexity-vs-fit Pareto front for discovered models.
- objective=:rss minimizes residual sum-of-squares proxy (residual_norm^2)
- objective=:r2 maximizes R^2 (internally minimized as -R^2)
"""
function compute_pareto_front(
    models::Vector{<:DiscoveredModel};
    objective::Symbol = :rss,
    r2_values::Union{Nothing, AbstractVector{<:Real}} = nothing,
)
    isempty(models) && return (
        indices = Int[],
        models = DiscoveredModel[],
        complexities = Int[],
        scores = Float64[],
        objective = objective,
    )

    scores, raw_scores = _model_objective_scores(models; objective = objective, r2_values = r2_values)
    complexities = [_model_complexity(m) for m in models]

    sorted_idx = sortperm(1:length(models); by = i -> (complexities[i], scores[i]))
    front_idx = Int[]
    best_score = Inf
    eps = 1e-12

    for i in sorted_idx
        if scores[i] < best_score - eps
            push!(front_idx, i)
            best_score = scores[i]
        end
    end

    return (
        indices = front_idx,
        models = models[front_idx],
        complexities = complexities[front_idx],
        scores = raw_scores[front_idx],
        objective = objective,
    )
end

function _fold_ranges(n::Int, k::Int)
    k_eff = min(max(k, 1), n)
    base = div(n, k_eff)
    remn = mod(n, k_eff)

    ranges = UnitRange{Int}[]
    start_idx = 1
    for fold in 1:k_eff
        fold_len = base + (fold <= remn ? 1 : 0)
        stop_idx = start_idx + fold_len - 1
        push!(ranges, start_idx:stop_idx)
        start_idx = stop_idx + 1
    end
    return ranges
end

"""
    kfold_cv_residual(model, blocks, score_fn; k=5)

Compute mean K-fold residual metric over observation blocks.
`score_fn` must accept `(model, train_blocks, val_blocks)` and return a scalar.
"""
function kfold_cv_residual(
    model::DiscoveredModel,
    blocks::Vector{ObservationTable},
    score_fn::Function;
    k::Int = 5,
)
    n = length(blocks)
    n > 0 || throw(ArgumentError("blocks cannot be empty."))

    folds = _fold_ranges(n, k)
    fold_scores = Float64[]

    for val_range in folds
        val_idx = collect(val_range)
        train_idx = [i for i in 1:n if !(i in val_range)]

        train_blocks = blocks[train_idx]
        val_blocks = blocks[val_idx]

        score = Float64(score_fn(model, train_blocks, val_blocks))
        isfinite(score) || throw(ArgumentError("score_fn returned a non-finite score."))
        push!(fold_scores, score)
    end

    return (
        mean_score = mean(fold_scores),
        fold_scores = fold_scores,
        k = length(fold_scores),
    )
end
