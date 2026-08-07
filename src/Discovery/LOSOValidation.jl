# src/Discovery/LOSOValidation.jl
using CSV
using DataFrames
using JSON3
using LinearAlgebra
using Printf
using Statistics

"""
    CampaignData

Container for one campaign/site used by LOSO cross-validation.

Fields:
- `site_name`: human-readable site identifier.
- `states`: sample matrix where rows are samples and columns are candidate features/targets.
- `time`: sample timestamps aligned with rows of `states`.
- `z`: vertical coordinates associated with observations.
- `auxiliary`: optional metadata for discovery/evaluation.
"""
struct CampaignData
    site_name::String
    states::Matrix{Float64}
    time::Vector{Float64}
    z::Vector{Float64}
    auxiliary::Dict{Symbol, Any}

    function CampaignData(
        site_name::String,
        states::AbstractMatrix{<:Real},
        time::AbstractVector{<:Real},
        z::AbstractVector{<:Real},
        auxiliary::Dict{Symbol, Any} = Dict{Symbol, Any}(),
    )
        ns = size(states, 1)
        length(time) == ns || throw(ArgumentError("time length must equal number of rows in states."))
        return new(
            site_name,
            Matrix{Float64}(states),
            Vector{Float64}(time),
            Vector{Float64}(z),
            auxiliary,
        )
    end
end

"""
    LOSOSplit

One Leave-One-Site-Out split with training campaigns and one held-out validation campaign.
"""
struct LOSOSplit
    train_campaigns::Vector{CampaignData}
    val_campaign::CampaignData
    val_site_name::String
end

"""
    LOSOResult

Per-site LOSO evaluation record.
"""
struct LOSOResult
    val_site_name::String
    discovered_model::String
    coefficients::Vector{Float64}
    train_rss::Float64
    val_rmse::Float64
    val_mae::Float64
    val_r2::Float64
    coverage_probability::Float64
end

"""
    LOSOSummary

Aggregate LOSO summary across all held-out sites.
"""
struct LOSOSummary
    results::Vector{LOSOResult}
    mean_validation_rmse::Float64
    cross_campaign_stability_score::Float64
end

function _escape_tex_text(str::AbstractString)
    escaped = String(str)
    escaped = replace(escaped, "\\" => "\\textbackslash{}")
    escaped = replace(escaped, "_" => "\\_")
    escaped = replace(escaped, "%" => "\\%")
    escaped = replace(escaped, "&" => "\\&")
    escaped = replace(escaped, "#" => "\\#")
    escaped = replace(escaped, "{" => "\\{")
    escaped = replace(escaped, "}" => "\\}")
    escaped = replace(escaped, string('$') => string('\\', '$'))
    return escaped
end

function _fmt4(x::Real)
    return isfinite(x) ? @sprintf("%.4f", Float64(x)) : "NA"
end

function _fmt3(x::Real)
    return isfinite(x) ? @sprintf("%.3f", Float64(x)) : "NA"
end

function _fmt_pct(x::Real)
    return isfinite(x) ? @sprintf("%.1f\\%%", 100.0 * Float64(x)) : "NA"
end

"""
    export_loso_table(summary::LOSOSummary, filepath::String)

Export a `LOSOSummary` object to a LaTeX tabular snippet matching `report.tex.mustache`.
"""
function export_loso_table(summary::LOSOSummary, filepath::String)
    mkpath(dirname(filepath))
    open(filepath, "w") do io
        dollar_tex = string('$')
        println(io, "\\begin{tabular}{lcccc}")
        println(io, "\\toprule")
        println(io, "Held-Out Site & Validation RMSE & Validation MAE & " * dollar_tex * "R^2" * dollar_tex * " & 95\\% Coverage \\\\")
        println(io, "\\midrule")

        for r in summary.results
            site = _escape_tex_text(r.val_site_name)
            println(io, "$(site) & $(_fmt4(r.val_rmse)) & $(_fmt4(r.val_mae)) & $(_fmt4(r.val_r2)) & $(_fmt_pct(r.coverage_probability)) \\\\")
        end

        println(io, "\\midrule")
        println(io, "\\multicolumn{5}{l}{\\textbf{Mean Validation RMSE}: $(_fmt4(summary.mean_validation_rmse)) \\quad \\textbf{Stability Score}: $(_fmt3(summary.cross_campaign_stability_score))} \\\\")
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
    end
    return filepath
end

function _site_slug(site_name::String)
    slug = lowercase(strip(site_name))
    slug = replace(slug, "-" => "_")
    slug = replace(slug, " " => "_")
    return slug
end

function _candidate_paths(base_dirs::Vector{String}, names::Vector{String})
    paths = String[]
    for d in base_dirs
        for n in names
            push!(paths, joinpath(d, n))
        end
    end
    return paths
end

function _first_existing_path(paths::Vector{String})
    for p in paths
        if isfile(p)
            return p
        end
    end
    return nothing
end

function _loso_to_f64(x)
    if x === missing || x === nothing
        return NaN
    elseif x isa Number
        return Float64(x)
    end
    v = tryparse(Float64, strip(String(x)))
    return v === nothing ? NaN : v
end

function _vector_from_csv(path::String)
    df = CSV.read(path, DataFrame)
    ncol(df) >= 1 || throw(ArgumentError("CSV file $(path) has no columns."))
    return _loso_to_f64.(df[!, 1])
end

function _load_optional_metadata(artifact_dir::String, slug::String)
    site_name_upper = uppercase(slug)
    site_name_mixed = replace(site_name_upper, "_" => "-")

    json_dirs = String[
        joinpath(artifact_dir, "json"),
        joinpath(artifact_dir, "campaign_exports", "json"),
        joinpath(artifact_dir, "pde_benchmark"),
        artifact_dir,
    ]

    candidate_jsons = _candidate_paths(json_dirs, String[
        "$(slug)_model_and_diagnostics.json",
        "$(slug).json",
        "$(site_name_upper).json",
        "$(site_name_mixed).json",
        "campaigns.json",
        "pde_benchmark.json",
        "benchmark_summary.json",
    ])

    metadata = Dict{Symbol, Any}()
    source = _first_existing_path(candidate_jsons)
    if source === nothing
        return metadata
    end

    # Some generated files can be empty when upstream runs fail; treat as no metadata.
    if filesize(source) == 0
        return metadata
    end

    try
        parsed = JSON3.read(read(source, String))
        metadata[:source] = source
        metadata[:payload] = parsed
    catch
        metadata[:source] = source
        metadata[:payload] = Dict{Symbol, Any}()
    end
    return metadata
end

function _load_profile_bundle(artifact_dir::String, slug::String)
    csv_dirs = String[
        joinpath(artifact_dir, "csv"),
        joinpath(artifact_dir, "campaign_exports", "csv"),
        artifact_dir,
    ]

    profile_path = _first_existing_path(_candidate_paths(csv_dirs, String[
        "$(slug)_profile.csv",
        "$(slug)_states.csv",
    ]))
    z_path = _first_existing_path(_candidate_paths(csv_dirs, String[
        "$(slug)_z.csv",
        "$(slug)_heights.csv",
    ]))
    t_path = _first_existing_path(_candidate_paths(csv_dirs, String[
        "$(slug)_time.csv",
        "$(slug)_timestamps.csv",
    ]))

    if profile_path === nothing || z_path === nothing || t_path === nothing
        return nothing
    end

    U = Matrix{Float64}(CSV.read(profile_path, DataFrame))
    z = _vector_from_csv(z_path)
    t = _vector_from_csv(t_path)

    size(U, 1) == length(z) || throw(DimensionMismatch("Profile rows must match z length for $(slug)."))
    size(U, 2) == length(t) || throw(DimensionMismatch("Profile columns must match time length for $(slug)."))

    states = permutedims(U) # store as samples x features for downstream LOSO model matrix usage
    aux = Dict{Symbol, Any}(
        :ingestion_mode => :profile_bundle,
        :profile_path => profile_path,
        :z_path => z_path,
        :time_path => t_path,
        :profile_matrix => U,
    )

    return (states = states, time = t, z = z, auxiliary = aux)
end

function _load_raw_campaign_table(artifact_dir::String, slug::String)
    site_name_upper = uppercase(slug)
    site_name_mixed = replace(site_name_upper, "_" => "-")

    csv_dirs = String[
        joinpath(artifact_dir, "csv"),
        joinpath(artifact_dir, "campaign_exports", "csv"),
        artifact_dir,
    ]
    raw_path = _first_existing_path(_candidate_paths(csv_dirs, String[
        "$(slug)_raw.csv",
        "$(slug).csv",
        "$(site_name_upper).csv",
        "$(site_name_mixed).csv",
        "$(slug)_profiles.csv",
    ]))

    if raw_path === nothing
        for (root, _, files) in walkdir(artifact_dir)
            for file in files
                if endswith(lowercase(file), ".csv") && contains(lowercase(file), slug)
                    raw_path = joinpath(root, file)
                    break
                end
            end
            raw_path === nothing || break
        end
    end

    raw_path === nothing && throw(ArgumentError("No campaign CSV artifact found for site slug $(slug)."))

    df = CSV.read(raw_path, DataFrame)
    nrow(df) > 0 || throw(ArgumentError("CSV artifact $(raw_path) is empty."))

    preferred_cols = Symbol[:zeta, :phi_obs, :phi_m, :phi_h, :L_obukhov, :ustar, :ws_lo, :ws_hi]
    present = [c for c in preferred_cols if hasproperty(df, c)]

    if length(present) < 2
        numeric_cols = names(df, Number)
        present = Symbol.(numeric_cols)
    end

    time_col = hasproperty(df, :time) ? :time : (hasproperty(df, :t) ? :t : nothing)
    z_col = hasproperty(df, :z) ? :z : (hasproperty(df, :height) ? :height : nothing)

    # For generic numeric-table ingestion, avoid treating time/height as model features.
    present = [c for c in present if c != time_col && c != z_col]

    length(present) >= 2 || throw(ArgumentError("Need at least two numeric state columns to build LOSO states for $(slug)."))

    M = hcat([_loso_to_f64.(df[!, c]) for c in present]...)

    time_vec = if time_col !== nothing
        _loso_to_f64.(df[!, time_col])
    else
        collect(1.0:nrow(df))
    end

    z_vec = if z_col !== nothing
        _loso_to_f64.(df[!, z_col])
    else
        collect(1.0:nrow(df))
    end

    finite_mask = vec(all(isfinite, M; dims = 2)) .& isfinite.(time_vec)
    if !any(finite_mask)
        throw(ArgumentError("All candidate rows are non-finite for site $(slug)."))
    end

    states = Matrix{Float64}(M[finite_mask, :])
    time = Vector{Float64}(time_vec[finite_mask])
    z = Vector{Float64}(z_vec[finite_mask])

    aux = Dict{Symbol, Any}(
        :ingestion_mode => :raw_table,
        :raw_csv_path => raw_path,
        :column_names => present,
        :dropped_rows => count(.!finite_mask),
    )

    return (states = states, time = time, z = z, auxiliary = aux)
end

"""
    load_campaign_data(artifact_dir, site_name) -> CampaignData

Load a single campaign artifact bundle into `CampaignData`.

Resolution strategy:
1. Try explicit profile bundle files (`*_profile.csv`, `*_z.csv`, `*_time.csv`).
2. Fall back to parsing `*_raw.csv` and constructing a finite feature matrix.
3. Attach optional JSON metadata from campaign/pde benchmark summaries into `auxiliary`.
"""
function load_campaign_data(artifact_dir::String, site_name::String)::CampaignData
    slug = _site_slug(site_name)

    loaded = _load_profile_bundle(artifact_dir, slug)
    if loaded === nothing
        loaded = _load_raw_campaign_table(artifact_dir, slug)
    end

    metadata = _load_optional_metadata(artifact_dir, slug)
    auxiliary = Dict{Symbol, Any}(loaded.auxiliary)
    auxiliary[:metadata] = get(metadata, :payload, Dict{Symbol, Any}())
    auxiliary[:metadata_source] = get(metadata, :source, "")

    payload = auxiliary[:metadata]
    try
        for (k, v) in pairs(payload)
            ks = Symbol(String(k))
            if !haskey(auxiliary, ks)
                auxiliary[ks] = v
            end
        end
    catch
        # Keep metadata nested-only when payload is not pair-iterable.
    end

    return CampaignData(site_name, loaded.states, loaded.time, loaded.z, auxiliary)
end

"""
    load_all_campaigns(artifacts_dir; sites=["SHEBA", "CASES-99", "FLOSS", "BLLAST"]) -> Vector{CampaignData}

Load all available campaign artifacts for the requested sites, warning and skipping
sites that cannot be loaded.
"""
function load_all_campaigns(
    artifacts_dir::String;
    sites::Vector{String} = ["SHEBA", "CASES-99", "FLOSS", "BLLAST"],
)::Vector{CampaignData}
    out = CampaignData[]
    for site in sites
        try
            push!(out, load_campaign_data(artifacts_dir, site))
        catch err
            @warn "Skipping site during artifact ingestion" site error = sprint(showerror, err)
        end
    end
    return out
end

"""
    create_loso_splits(campaigns)

Generate all Leave-One-Site-Out partitions from campaign data.
"""
function create_loso_splits(campaigns::Vector{CampaignData})::Vector{LOSOSplit}
    n = length(campaigns)
    n >= 2 || throw(ArgumentError("LOSO requires at least two campaigns."))

    splits = Vector{LOSOSplit}(undef, n)
    for i in 1:n
        val = campaigns[i]
        train = CampaignData[campaigns[j] for j in 1:n if j != i]
        splits[i] = LOSOSplit(train, val, val.site_name)
    end
    return splits
end

function _r2_score(pred::Vector{Float64}, truth::Vector{Float64})
    length(pred) == length(truth) || throw(DimensionMismatch("pred and truth must match in length."))
    μ = mean(truth)
    rss = sum(abs2, pred .- truth)
    tss = sum(abs2, truth .- μ)
    return tss <= 1e-12 ? 0.0 : 1.0 - rss / tss
end

"""
    evaluate_out_of_sample(model_coeffs, library_matrix_val, target_val)

Evaluate out-of-sample performance on validation data.

Returns a named tuple containing:
- `rmse`
- `mae`
- `r2`
- `coverage_probability` (empirical coverage under a Gaussian residual 95% interval)
"""
function evaluate_out_of_sample(
    model_coeffs::AbstractVector{<:Real},
    library_matrix_val::AbstractMatrix{<:Real},
    target_val::AbstractVector{<:Real},
)
    X = Matrix{Float64}(library_matrix_val)
    y = Vector{Float64}(target_val)
    β = Vector{Float64}(model_coeffs)

    size(X, 1) == length(y) || throw(DimensionMismatch("library_matrix_val rows must match target length."))
    size(X, 2) == length(β) || throw(DimensionMismatch("model_coeffs length must match library_matrix_val columns."))

    ŷ = X * β
    resid = ŷ .- y

    rmse = sqrt(mean(abs2, resid))
    mae = mean(abs.(resid))
    r2 = _r2_score(ŷ, y)

    σ = std(resid)
    if !isfinite(σ) || σ <= 0.0
        coverage = mean(abs.(resid) .<= 1e-10)
    else
        lower = ŷ .- 1.96 * σ
        upper = ŷ .+ 1.96 * σ
        coverage = mean((y .>= lower) .& (y .<= upper))
    end

    return (rmse = rmse, mae = mae, r2 = r2, coverage_probability = coverage)
end

function _campaign_design_matrix(
    campaign::CampaignData;
    feature_cols::Union{Nothing, AbstractVector{Int}} = nothing,
    target_col::Union{Nothing, Int} = nothing,
)
    ncols = size(campaign.states, 2)
    ncols >= 2 || throw(ArgumentError("campaign.states must contain at least one feature column and one target column."))

    tgt = isnothing(target_col) ? ncols : target_col
    1 <= tgt <= ncols || throw(ArgumentError("target_col out of bounds."))

    feats = isnothing(feature_cols) ? [j for j in 1:ncols if j != tgt] : collect(feature_cols)
    isempty(feats) && throw(ArgumentError("feature_cols cannot be empty."))
    all(1 <= j <= ncols for j in feats) || throw(ArgumentError("feature_cols contains out-of-bounds index."))
    any(j -> j == tgt, feats) && throw(ArgumentError("feature_cols cannot include target_col."))

    X = campaign.states[:, feats]
    y = campaign.states[:, tgt]
    return Matrix{Float64}(X), Vector{Float64}(y)
end

function _stack_training_data(
    campaigns::Vector{CampaignData};
    feature_cols::Union{Nothing, AbstractVector{Int}} = nothing,
    target_col::Union{Nothing, Int} = nothing,
)
    Xs = Matrix{Float64}[]
    ys = Vector{Float64}[]
    for c in campaigns
        X, y = _campaign_design_matrix(c; feature_cols = feature_cols, target_col = target_col)
        push!(Xs, X)
        push!(ys, y)
    end
    return vcat(Xs...), vcat(ys...)
end

function _extract_coefficients(output)
    if output isa NamedTuple
        if haskey(output, :coefficients)
            return Vector{Float64}(output.coefficients)
        end
    end

    if output isa Tuple && !isempty(output)
        first_item = output[1]
        if first_item isa AbstractVector{<:Real}
            return Vector{Float64}(first_item)
        end
    end

    if output isa AbstractVector{<:Real}
        return Vector{Float64}(output)
    end

    throw(ArgumentError("discovery_fn must return coefficients as a vector, a tuple with vector first element, or a named tuple with :coefficients."))
end

function _extract_model_string(output)
    if output isa NamedTuple && haskey(output, :discovered_model)
        return string(output.discovered_model)
    elseif output isa NamedTuple && haskey(output, :model)
        return string(output.model)
    end
    return string(typeof(output))
end

function _extract_train_rss(output, Xtrain::Matrix{Float64}, ytrain::Vector{Float64}, coeffs::Vector{Float64})
    if output isa NamedTuple && haskey(output, :train_rss)
        return float(output.train_rss)
    end
    resid = Xtrain * coeffs .- ytrain
    return sum(abs2, resid)
end

function _stability_score(results::Vector{LOSOResult})
    n = length(results)
    n <= 1 && return 1.0

    norms = [norm(r.coefficients) for r in results]
    sims = Float64[]

    for i in 1:(n - 1)
        for j in (i + 1):n
            ci = results[i].coefficients
            cj = results[j].coefficients
            length(ci) == length(cj) || continue

            ni = norms[i]
            nj = norms[j]
            if ni <= 1e-12 || nj <= 1e-12
                push!(sims, 0.0)
            else
                cos_sim = dot(ci, cj) / (ni * nj)
                push!(sims, clamp((1.0 + cos_sim) / 2.0, 0.0, 1.0))
            end
        end
    end

    return isempty(sims) ? 0.0 : mean(sims)
end

"""
    run_loso_cross_validation(campaigns, discovery_fn; kwargs...)

Run Leave-One-Site-Out cross-validation.

`discovery_fn` is called as:
`discovery_fn(train_campaigns, X_train, y_train; kwargs...)`
and must return coefficients via one of:
- `NamedTuple` with `:coefficients`
- tuple with coefficient vector as first element
- coefficient vector directly

Keyword arguments:
- `feature_cols`: feature column indices from `CampaignData.states` (default: all except target)
- `target_col`: target column index from `CampaignData.states` (default: last column)

Returns a `LOSOSummary`.
"""
function run_loso_cross_validation(
    campaigns::Vector{CampaignData},
    discovery_fn::Function;
    feature_cols::Union{Nothing, AbstractVector{Int}} = nothing,
    target_col::Union{Nothing, Int} = nothing,
    kwargs...
)
    splits = create_loso_splits(campaigns)
    results = LOSOResult[]

    for split in splits
        X_train, y_train = _stack_training_data(
            split.train_campaigns;
            feature_cols = feature_cols,
            target_col = target_col,
        )
        X_val, y_val = _campaign_design_matrix(
            split.val_campaign;
            feature_cols = feature_cols,
            target_col = target_col,
        )

        raw = discovery_fn(split.train_campaigns, X_train, y_train; kwargs...)
        coeffs = _extract_coefficients(raw)
        size(X_val, 2) == length(coeffs) || throw(DimensionMismatch("Validation design matrix columns must match discovered coefficient length."))

        metrics = evaluate_out_of_sample(coeffs, X_val, y_val)
        train_rss = _extract_train_rss(raw, X_train, y_train, coeffs)
        model_str = _extract_model_string(raw)

        push!(results, LOSOResult(
            split.val_site_name,
            model_str,
            coeffs,
            train_rss,
            metrics.rmse,
            metrics.mae,
            metrics.r2,
            metrics.coverage_probability,
        ))
    end

    mean_rmse = mean(r.val_rmse for r in results)
    stability = _stability_score(results)
    return LOSOSummary(results, mean_rmse, stability)
end

"""
    run_artifact_loso(artifacts_dir, discovery_fn; kwargs...) -> LOSOSummary

Convenience wrapper that ingests campaign artifacts and runs LOSO cross-validation
in one call.

Accepted keyword arguments include:
- `sites`
- any `run_loso_cross_validation` keyword arguments (e.g. `feature_cols`, `target_col`)
"""
function run_artifact_loso(
    artifacts_dir::String,
    discovery_fn::Function;
    sites::Vector{String} = ["SHEBA", "CASES-99", "FLOSS", "BLLAST"],
    kwargs...
)::LOSOSummary
    campaigns = load_all_campaigns(artifacts_dir; sites = sites)
    length(campaigns) >= 2 || throw(ArgumentError("Need at least two ingested campaigns to run LOSO."))
    return run_loso_cross_validation(campaigns, discovery_fn; kwargs...)
end
