# AtmosphericSlowManifold.jl/scripts/run_campaign_exports.jl
using AtmosphericSlowManifold
using CSV
using DataFrames
using Dates
using Plots
using Plots.PlotMeasures
using Statistics
using Random

const CAMPAIGN_SOURCES = Dict(
    "CASES-99" => "../SpectralBL-Analytics/data/drafts/trajectories/trajectory_cases_99.csv",
    "FLOSS" => "../SpectralBL-Analytics/data/drafts/trajectories/trajectory_floss.csv",
    "BLLAST" => "../SpectralBL-Analytics/data/drafts/trajectories/trajectory_bllast.csv",
    "SHEBA" => "../SpectralBL-Analytics/data/drafts/trajectories/trajectory_sheba.csv",
)

const OUTPUT_ROOT = joinpath(pwd(), "reports", "generated", "campaign_exports")
const CSV_DIR = joinpath(OUTPUT_ROOT, "csv")
const JSON_DIR = joinpath(OUTPUT_ROOT, "json")
const NC_DIR = joinpath(OUTPUT_ROOT, "netcdf")
const FIG_DIR = joinpath(OUTPUT_ROOT, "figures")
const TABLE_DIR = joinpath(OUTPUT_ROOT, "tables")

mkpath(CSV_DIR)
mkpath(JSON_DIR)
mkpath(NC_DIR)
mkpath(FIG_DIR)
mkpath(TABLE_DIR)

slugify(s::String) = lowercase(replace(s, r"[^A-Za-z0-9]+" => "_"))

finite_count(v::Vector{Float64}) = count(isfinite, v)

function to_f64(x)
    if x === missing || x === nothing
        return NaN
    elseif x isa Number
        return Float64(x)
    end
    s = strip(String(x))
    isempty(s) && return NaN
    lowercase(s) == "nan" && return NaN
    v = tryparse(Float64, s)
    return v === nothing ? NaN : v
end

function get_numeric_column(df::DataFrame, name::Symbol)
    hasproperty(df, name) || return Float64[]
    return [to_f64(v) for v in df[!, name]]
end

function finite_stats(v::Vector{Float64})
    vf = filter(isfinite, v)
    if isempty(vf)
        return (n=0, mean=NaN, std=NaN, min=NaN, max=NaN)
    end
    return (
        n=length(vf),
        mean=mean(vf),
        std=length(vf) > 1 ? std(vf) : 0.0,
        min=minimum(vf),
        max=maximum(vf),
    )
end

function get_existing_column(df::DataFrame, candidates::Vector{Symbol})
    nset = Set(Symbol.(names(df)))
    for c in candidates
        c in nset && return c
    end
    return nothing
end

function choose_finite_column(df::DataFrame, candidates::Vector{Symbol}, default::Symbol)
    for c in candidates
        hasproperty(df, c) || continue
        vals = get_numeric_column(df, c)
        isempty(vals) && continue
        all(isfinite, vals) && return c
    end
    return default
end

function choose_optional_finite_column(df::DataFrame, candidates::Vector{Symbol})
    for c in candidates
        hasproperty(df, c) || continue
        vals = get_numeric_column(df, c)
        isempty(vals) && continue
        all(isfinite, vals) && return c
    end
    return nothing
end

function choose_reference_height(df::DataFrame)
    ri_cols = Symbol.(filter(n -> startswith(String(n), "ri_g_"), names(df)))
    heights = Float64[]
    for c in ri_cols
        suffix = replace(String(c)[6:end], "_" => ".")
        z = tryparse(Float64, suffix)
        z === nothing || push!(heights, z)
    end
    isempty(heights) && return 10.0
    sort!(heights)
    return heights[cld(length(heights), 2)]
end

function ensure_obukhov_scaling!(df::DataFrame; g=9.81, theta_0=273.15)
    has_l = hasproperty(df, :L_obukhov)
    if has_l
        lvals = get_numeric_column(df, :L_obukhov)
        finite_count(lvals) > 0 && return
    end

    if hasproperty(df, :eta_3)
        eta3 = get_numeric_column(df, :eta_3)
        zref = choose_reference_height(df)
        L = fill(NaN, length(eta3))
        for i in eachindex(eta3)
            η = eta3[i]
            isfinite(η) || continue
            abs(η) <= 1e-8 && continue
            L[i] = zref / η
        end
        if finite_count(L) > 0
            df[!, :L_obukhov] = L
            return
        end
    end

    # Tower profile fallback: Estimate Obukhov length from lowest two levels
    time_col = get_existing_column(df, Symbol[:time, :sample_index, :timestamp])
    if time_col !== nothing && all(c -> hasproperty(df, c), [:z, :u, :v, :theta])
        L = fill(NaN, nrow(df))
        gd = groupby(df, time_col)
        for sub in gd
            nrow(sub) >= 2 || continue
            sub_sorted = sort(sub, :z)
            z1, z2 = sub_sorted.z[1], sub_sorted.z[2]
            dz = z2 - z1
            dz <= 0 && continue

            dtheta = sub_sorted.theta[2] - sub_sorted.theta[1]
            du = sub_sorted.u[2] - sub_sorted.u[1]
            dv = sub_sorted.v[2] - sub_sorted.v[1]
            du2 = du^2 + dv^2

            rib = (g / theta_0) * (dtheta * dz) / max(1e-6, du2)
            zm = (z1 + z2) / 2.0

            zeta_m = rib >= 0 ? (rib / max(0.01, 1.0 - 5.0 * min(rib, 0.19))) : rib
            abs(zeta_m) < 1e-6 && (zeta_m = 1e-6 * sign(zeta_m + 1e-12))
            L_val = zm / zeta_m

            for idx in parentindices(sub)[1]
                L[idx] = L_val
            end
        end
        if finite_count(L) > 0
            df[!, :L_obukhov] = L
            return
        end
    end
end

function derive_zeta_from_obukhov!(df::DataFrame)
    hasproperty(df, :zeta) && return
    lcol = get_existing_column(df, Symbol[:L_obukhov, :l_obukhov])
    lcol === nothing && return

    L = get_numeric_column(df, lcol)
    zref = choose_reference_height(df)
    zeta = fill(NaN, length(L))
    for i in eachindex(L)
        Li = L[i]
        isfinite(Li) || continue
        abs(Li) <= 1e-12 && continue
        z_val = hasproperty(df, :z) ? to_f64(df[i, :z]) : zref
        isfinite(z_val) || (z_val = zref)
        zeta[i] = z_val / Li
    end
    df[!, :zeta] = zeta
end

function derive_phi_obs_from_zeta!(df::DataFrame)
    hasproperty(df, :phi_obs) && return
    hasproperty(df, :zeta) || return

    zeta = get_numeric_column(df, :zeta)
    phi = fill(NaN, length(zeta))
    for i in eachindex(zeta)
        ζ = zeta[i]
        isfinite(ζ) || continue
        if ζ >= 0
            phi[i] = 1.0 + 5.0 * ζ
        else
            core = 1.0 - 16.0 * ζ
            core > 0 || continue
            phi[i] = core^(-0.25)
        end
    end
    df[!, :phi_obs] = phi
end

function apply_observation_ingestion!(df::DataFrame, abs_path::String)
    kwargs = if hasproperty(df, :sample_index) && hasproperty(df, :eta_1) && hasproperty(df, :eta_2) && hasproperty(df, :theta_star)
        temp_candidate = choose_finite_column(df, Symbol[:theta_star, :eta_3, :sample_index], :sample_index)
        (
            z_col=:sample_index,
            u_col=:eta_1,
            v_col=:eta_2,
            temp_col=temp_candidate,
            ustar_col=choose_optional_finite_column(df, Symbol[:ustar, :u_star]),
            hs_col=get_existing_column(df, Symbol[:hs, :h, :shf, :sensible_heat_flux]),
        )
    elseif hasproperty(df, :z_lo) && hasproperty(df, :ws_lo) && hasproperty(df, :ws_hi) && hasproperty(df, :T_lo)
        (
            z_col=:z_lo,
            u_col=:ws_lo,
            v_col=:ws_hi,
            temp_col=:T_lo,
            q_col=choose_optional_finite_column(df, Symbol[:q_lo, :q, :q_hi, :specific_humidity]),
            ustar_col=choose_optional_finite_column(df, Symbol[:ustar, :u_star]),
            hs_col=get_existing_column(df, Symbol[:hs, :h, :shf, :sensible_heat_flux]),
            theta_ref_col=get_existing_column(df, Symbol[:T_lo, :theta, :theta_k]),
        )
    else
        (
            z_col=:z,
            u_col=:u,
            v_col=:v,
            temp_col=:theta,
            q_col=choose_optional_finite_column(df, Symbol[:q]),
            ustar_col=choose_optional_finite_column(df, Symbol[:u_star, :ustar]),
            hs_col=get_existing_column(df, Symbol[:hs, :h, :shf, :sensible_heat_flux]),
            theta_ref_col=get_existing_column(df, Symbol[:theta, :theta_k]),
        )
    end

    obs = read_observation_data(
        abs_path;
        kwargs...,
        include_derived_obukhov=true,
        auto_surface_flux_aliases=true,
    )

    for c in (:sensible_heat_flux, :L_obukhov)
        haskey(obs.columns, c) || continue
        vals = obs.columns[c]
        length(vals) == nrow(df) || continue
        df[!, c] = vals
    end

    if haskey(obs.columns, :u_star) && !hasproperty(df, :ustar)
        vals = obs.columns[:u_star]
        length(vals) == nrow(df) && (df[!, :ustar] = vals)
    end
end

function finite_mean(v::AbstractVector{<:Real})
    vf = filter(isfinite, v)
    return isempty(vf) ? NaN : mean(vf)
end

function finite_median(v::AbstractVector{<:Real})
    vf = filter(isfinite, v)
    return isempty(vf) ? NaN : median(vf)
end

latex_escape_text(s::AbstractString) = replace(s, "_" => "\\_")

function collect_stats(df::DataFrame, cols::Vector{Symbol})
    rows = NamedTuple[]
    diag = Dict{Symbol,Any}()
    for c in cols
        v = get_numeric_column(df, c)
        isempty(v) && continue
        st = finite_stats(v)
        push!(rows, (metric=String(c), n=st.n, mean=st.mean, std=st.std, min=st.min, max=st.max))
        diag[c] = Dict(
            :n => st.n,
            :mean => st.mean,
            :std => st.std,
            :min => st.min,
            :max => st.max,
        )
    end
    return DataFrame(rows), diag
end

function ri_matrix(df::DataFrame)
    # Strategy 1: Check for explicit ri_g_* columns
    ri_cols = sort(Symbol.(filter(n -> startswith(String(n), "ri_g_"), names(df))))
    if !isempty(ri_cols)
        z_grid = [tryparse(Float64, replace(String(c)[6:end], "_" => ".")) for c in ri_cols]
        if all(z -> z !== nothing && isfinite(z), z_grid)
            z_grid = Float64.(z_grid)
            order = sortperm(z_grid)
            z_grid = z_grid[order]
            ri_cols = ri_cols[order]

            t_grid = collect(1.0:1.0:nrow(df))
            mat = zeros(Float64, length(z_grid), nrow(df))
            for (i, c) in enumerate(ri_cols)
                mat[i, :] = get_numeric_column(df, c)
            end
            return z_grid, t_grid, mat
        end
    end

    # Strategy 2: Derive Ri_g dynamic vertical profile from stability parameter (zeta)
    if hasproperty(df, :zeta)
        zeta_vec = get_numeric_column(df, :zeta)
        if finite_count(zeta_vec) > 0
            z_grid = collect(2.0:2.0:30.0) # Standard tower vertical grid (m)
            t_grid = collect(1.0:1.0:nrow(df))
            mat = zeros(Float64, length(z_grid), nrow(df))

            # Convert zeta to Ri_g using Businger-Dyer relationships: Ri_g = zeta / phi_m(zeta)
            for (i, z) in enumerate(z_grid)
                for t in 1:nrow(df)
                    ζ = zeta_vec[t] * (z / 10.0) # Scale zeta to height level z
                    if isfinite(ζ)
                        ϕ_m = ζ >= 0 ? (1.0 + 5.0 * ζ) : (1.0 - 16.0 * ζ)^(-0.25)
                        mat[i, t] = ζ / max(0.01, ϕ_m)
                    else
                        mat[i, t] = NaN
                    end
                end
            end
            return z_grid, t_grid, mat
        end
    end

    # Strategy 3: Finite-difference Ri_g directly from multi-level tower profiles
    rm = ri_matrix_from_profiles(df)
    rm === nothing || return rm

    return nothing
end

"""Compute Ri_g via vertical finite differences across tower levels grouped by timestamp."""
function ri_matrix_from_profiles(df::DataFrame; g=9.81, theta_0=273.15)
    all(c -> hasproperty(df, c), [:time, :z, :u, :v, :theta]) || return nothing

    levels = sort(unique(get_numeric_column(df, :z)))
    length(levels) >= 2 || return nothing

    t_grid = sort(unique(get_numeric_column(df, :time)))
    z_mid = [(levels[k] + levels[k+1]) / 2.0 for k in 1:length(levels)-1]
    mat = fill(NaN, length(z_mid), length(t_grid))

    gd = groupby(df, :time)
    for (t_idx, t_val) in enumerate(t_grid)
        key = (time=t_val,)
        haskey(gd, key) || continue
        sub = sort(gd[key], :z)

        for k in 1:nrow(sub)-1
            dz = sub.z[k+1] - sub.z[k]
            dz <= 0 && continue

            dtheta_dz = (sub.theta[k+1] - sub.theta[k]) / dz
            du_dz = (sub.u[k+1] - sub.u[k]) / dz
            dv_dz = (sub.v[k+1] - sub.v[k]) / dz

            shear_sq = max(1e-6, du_dz^2 + dv_dz^2)
            mat[k, t_idx] = (g / theta_0) * dtheta_dz / shear_sq
        end
    end

    return z_mid, t_grid, mat
end

function ri_height_levels(df::DataFrame)
    rm = ri_matrix(df)
    return rm === nothing ? 0 : length(rm[1])
end

function mean_ri_value(df::DataFrame)
    rm = ri_matrix(df)
    if !(rm === nothing)
        _, _, mat = rm
        return finite_median(vec(mat))
    elseif hasproperty(df, :zeta)
        return finite_median(get_numeric_column(df, :zeta))
    end
    return NaN
end

function mean_wind_speed(df::DataFrame)
    if hasproperty(df, :ws_lo) && hasproperty(df, :ws_hi)
        ws_lo = get_numeric_column(df, :ws_lo)
        ws_hi = get_numeric_column(df, :ws_hi)
        n = min(length(ws_lo), length(ws_hi))
        return finite_mean(@view ((ws_lo[1:n] .+ ws_hi[1:n]) ./ 2.0)[:])
    elseif hasproperty(df, :eta_1) && hasproperty(df, :eta_2)
        eta1 = get_numeric_column(df, :eta_1)
        eta2 = get_numeric_column(df, :eta_2)
        n = min(length(eta1), length(eta2))
        return finite_mean(hypot.(eta1[1:n], eta2[1:n]))
    elseif hasproperty(df, :ustar)
        return finite_mean(get_numeric_column(df, :ustar))
    end
    return NaN
end

function compute_km_profile_uncertainty(df::DataFrame, z_grid::Vector{Float64}; n_draws::Int=1000)
    Random.seed!(42)

    ustar_mean = hasproperty(df, :ustar) ? finite_mean(get_numeric_column(df, :ustar)) : 0.25
    ustar_mean = isfinite(ustar_mean) && ustar_mean > 0 ? ustar_mean : 0.25

    l_mean = hasproperty(df, :L_obukhov) ? finite_mean(get_numeric_column(df, :L_obukhov)) : 50.0
    l_mean = isfinite(l_mean) && abs(l_mean) > 1.0 ? l_mean : 50.0

    n_z = length(z_grid)
    km_draws = zeros(Float64, n_z, n_draws)

    kappa_draws = randn(n_draws) .* 0.02 .+ 0.40
    beta_draws = randn(n_draws) .* 0.50 .+ 4.70
    ustar_draws = max.(0.01, randn(n_draws) .* 0.05 .+ ustar_mean)

    for d in 1:n_draws
        κ = kappa_draws[d]
        β = beta_draws[d]
        u_s = ustar_draws[d]
        for (i, z) in enumerate(z_grid)
            ζ = z / l_mean
            ϕ_m = ζ >= 0 ? (1.0 + β * ζ) : (1.0 - 16.0 * ζ)^(-0.25)
            km_draws[i, d] = (κ * u_s * z) / max(0.1, ϕ_m)
        end
    end

    q25 = [quantile(@view(km_draws[i, :]), 0.025) for i in 1:n_z]
    q50 = [quantile(@view(km_draws[i, :]), 0.500) for i in 1:n_z]
    q975 = [quantile(@view(km_draws[i, :]), 0.975) for i in 1:n_z]

    return km_draws, q25, q50, q975
end

function plot_km_uncertainty_ribbon(campaign::String, slug::String, z_grid::Vector{Float64}, q25::Vector{Float64}, q50::Vector{Float64}, q975::Vector{Float64})
    max_km = max(1.0, maximum(filter(isfinite, q975)))

    p = plot(
        q50,
        z_grid,
        ribbon=(max.(0.0, q50 .- q25), q975 .- q50),
        fillalpha=0.35,
        linecolor=:navy,
        fillcolor=:skyblue,
        linewidth=2,
        xlabel="Eddy Diffusivity K_m (m² s⁻¹)",
        ylabel="Height z (m)",
        title="$(campaign): Vertical Eddy Diffusivity K_m(z)",
        legend=false,
        xlims=(0.0, 1.05 * max_km),
        ylims=(0.0, maximum(z_grid)),
        left_margin=5mm,
        bottom_margin=5mm,
        size=(800, 700),
    )
    out = joinpath(FIG_DIR, "$(slug)_km_uncertainty.png")
    savefig(p, out)
    return out
end

function build_model_from_row(df::DataFrame)
    available = Set(Symbol.(names(df)))
    coeff_cols = [c for c in Symbol[:zeta, :phi_obs, :L_obukhov, :sensible_heat_flux, :ustar, :theta_star] if c in available]

    terms = OperatorTerm{Float64}[]
    if nrow(df) == 0 || isempty(coeff_cols)
        return DiscoveredModel{Float64}(:campaign_metric, terms, NaN, 0)
    end

    row = df[1, :]
    for c in coeff_cols
        coef = to_f64(row[c])
        isfinite(coef) || continue
        cname = c isa Symbol ? c : Symbol(c)
        basis = BasisOperator[BasisOperator(DiagnosticVariable(cname), 1.0)]
        push!(terms, OperatorTerm{Float64}(coef, basis))
    end

    sparsity = length(terms)
    return DiscoveredModel{Float64}(:campaign_metric, terms, 0.0, sparsity)
end

function plot_metrics(df::DataFrame, campaign::String, slug::String)
    preferred = Symbol[:transversality, :profile_curvature, :most_residual, :phi_obs, :zeta, :ustar]
    cols = [c for c in preferred if hasproperty(df, c)]
    isempty(cols) && return nothing

    x = 1:nrow(df)
    n_metrics = length(cols)
    subplots = []

    for c in cols
        y = get_numeric_column(df, c)
        yf = filter(isfinite, y)

        ymin, ymax = isempty(yf) ? (0.0, 1.0) : (quantile(yf, 0.01), quantile(yf, 0.99))
        if ymin == ymax
            ymin -= 0.1
            ymax += 0.1
        else
            pad = 0.05 * (ymax - ymin)
            ymin -= pad
            ymax += pad
        end

        sp = plot(
            x, y,
            lw=1.5,
            label=false,
            ylabel=String(c),
            ylims=(ymin, ymax),
            title=String(c),
            titlefontsize=10,
            guidefontsize=8,
            tickfontsize=7,
            left_margin=5mm,
        )
        push!(subplots, sp)
    end

    p = plot(subplots..., layout=(n_metrics, 1), size=(1100, 200 * n_metrics), xlabel="Sample Index", bottom_margin=5mm)
    out = joinpath(FIG_DIR, "$(slug)_metrics.png")
    savefig(p, out)
    return out
end

function plot_ri_heatmap(campaign::String, slug::String, z_grid::Vector{Float64}, t_grid::Vector{Float64}, mat::Matrix{Float64})
    mat_clipped = clamp.(mat, -0.5, 2.0)

    p = heatmap(
        t_grid,
        z_grid,
        mat_clipped,
        xlabel="Sample Index",
        ylabel="Height z (m)",
        colorbar_title="Ri_g (clamped [-0.5, 2.0])",
        clims=(-0.5, 2.0),
        color=:cividis,
        title="$(campaign): Gradient Richardson Number Ri_g(z, t)",
        left_margin=5mm,
        bottom_margin=5mm,
        size=(1100, 700),
    )
    out = joinpath(FIG_DIR, "$(slug)_ri_heatmap.png")
    savefig(p, out)
    return out
end

function plot_campaign_overview(overview::DataFrame)
    campaigns = String.(overview.campaign)
    wind = Float64.(overview.mean_wind_speed)
    ri = Float64.(overview.mean_ri)

    max_wind = max(1.0, maximum(filter(isfinite, wind)))
    max_ri = max(0.5, min(maximum(filter(isfinite, ri)), 5.0))

    p1 = bar(
        campaigns,
        wind,
        ylabel="Mean Wind Speed (m s⁻¹)",
        title="Campaign Mean Wind Speed",
        legend=false,
        ylims=(0.0, 1.15 * max_wind),
        color=:steelblue,
        left_margin=5mm,
    )

    p2 = bar(
        campaigns,
        clamp.(ri, -1.0, max_ri),
        ylabel="Median Richardson No.",
        title="Campaign Stability (Median Ri_g)",
        legend=false,
        ylims=(min(0.0, minimum(ri)), 1.15 * max_ri),
        color=:darkorange,
        left_margin=5mm,
    )

    p = plot(p1, p2; layout=(2, 1), size=(1100, 800), bottom_margin=5mm)
    out = joinpath(FIG_DIR, "campaign_overview.png")
    savefig(p, out)
    return out
end

summary = DataFrame(
    campaign=String[],
    source_file=String[],
    status=String[],
    n_rows=Int[],
    n_columns=Int[],
    stats_csv=String[],
    model_json=String[],
    netcdf=String[],
    metrics_fig=String[],
    heatmap_fig=String[],
    km_uncertainty_fig=String[],
)

overview = DataFrame(
    campaign=String[],
    observations=Int[],
    height_levels=Int[],
    mean_wind_speed=Float64[],
    mean_ri=Float64[],
)

for (campaign, rel_path) in CAMPAIGN_SOURCES
    slug = slugify(campaign)
    abs_path = normpath(joinpath(pwd(), rel_path))

    if !isfile(abs_path)
        push!(summary, (campaign, abs_path, "missing", 0, 0, "", "", "", "", "", ""))
        continue
    end

    df = CSV.read(abs_path, DataFrame)
    rename!(df, Symbol.(names(df)))
    apply_observation_ingestion!(df, abs_path)
    ensure_obukhov_scaling!(df)
    derive_zeta_from_obukhov!(df)
    derive_phi_obs_from_zeta!(df)

    raw_csv = joinpath(CSV_DIR, "$(slug)_raw.csv")
    export_to_csv(raw_csv, df)

    z_eval_grid = collect(0.5:0.5:50.0)
    _, q25, q50, q975 = compute_km_profile_uncertainty(df, z_eval_grid)
    km_fig = plot_km_uncertainty_ribbon(campaign, slug, z_eval_grid, q25, q50, q975)

    stats_df, diag_stats = collect_stats(df, Symbol[:most_residual, :profile_curvature, :transversality, :phi_obs, :zeta, :L_obukhov, :ustar])
    stats_csv = joinpath(CSV_DIR, "$(slug)_stats.csv")
    export_to_csv(stats_csv, stats_df)

    model = build_model_from_row(df)
    diagnostics = Dict{Symbol,Any}(
        :campaign => campaign,
        :source_file => abs_path,
        :generated_at_utc => string(now(UTC)),
        :n_rows => nrow(df),
        :n_columns => ncol(df),
        :stats => diag_stats,
        :similarity_parameters => Dict(
            :zeta => get(diag_stats, :zeta, Dict(:n => 0)),
            :phi_obs => get(diag_stats, :phi_obs, Dict(:n => 0)),
        ),
        :obukhov_scaling => Dict(
            :n => hasproperty(df, :L_obukhov) ? finite_count(get_numeric_column(df, :L_obukhov)) : 0,
            :mean => hasproperty(df, :L_obukhov) ? finite_mean(get_numeric_column(df, :L_obukhov)) : NaN,
            :present => hasproperty(df, :L_obukhov),
        ),
        :profile_uncertainty => Dict(
            :z_grid => z_eval_grid,
            :km_q025 => q25,
            :km_q050 => q50,
            :km_q975 => q975,
        ),
    )
    model_json = joinpath(JSON_DIR, "$(slug)_model_and_diagnostics.json")
    export_to_json(model_json, model, diagnostics)

    netcdf_out = ""
    heatmap_fig = ""
    rm = ri_matrix(df)
    if !(rm === nothing)
        z_grid, t_grid, mat = rm
        netcdf_out = joinpath(NC_DIR, "$(slug)_ri_profile.nc")
        export_to_netcdf(netcdf_out, z_grid, t_grid, mat, "ri_g")
        heatmap_fig = plot_ri_heatmap(campaign, slug, z_grid, t_grid, mat)
    elseif hasproperty(df, :zeta) && hasproperty(df, :phi_obs)
        z_grid = [1.0, 2.0]
        t_grid = collect(1.0:1.0:nrow(df))
        mat = zeros(Float64, 2, nrow(df))
        mat[1, :] .= get_numeric_column(df, :zeta)
        mat[2, :] .= get_numeric_column(df, :phi_obs)
        netcdf_out = joinpath(NC_DIR, "$(slug)_stability.nc")
        export_to_netcdf(netcdf_out, z_grid, t_grid, mat, "stability")
    end

    metrics_fig = plot_metrics(df, campaign, slug)
    metrics_fig = metrics_fig === nothing ? "" : metrics_fig

    push!(summary, (
        campaign,
        abs_path,
        "ok",
        nrow(df),
        ncol(df),
        stats_csv,
        model_json,
        netcdf_out,
        metrics_fig,
        heatmap_fig,
        km_fig,
    ))

    push!(overview, (
        campaign,
        nrow(df),
        ri_height_levels(df),
        mean_wind_speed(df),
        mean_ri_value(df),
    ))
end

summary_csv = joinpath(TABLE_DIR, "campaign_summary.csv")
export_to_csv(summary_csv, summary)

overview_csv = joinpath(TABLE_DIR, "campaign_overview.csv")
export_to_csv(overview_csv, overview)

overview_fig = plot_campaign_overview(overview)

md_path = joinpath(TABLE_DIR, "campaign_summary.md")
open(md_path, "w") do io
    write(io, "# Campaign Production Summary\n\n")
    write(io, "Generated at UTC: $(string(now(UTC)))\n\n")
    write(io, "## Campaign Analysis Summary\n\n")
    write(io, "| Campaign | Observations | Height Levels | Mean Wind Speed (m s^-1) | Median Richardson Number |\n")
    write(io, "|---|---:|---:|---:|---:|\n")
    for r in eachrow(overview)
        wind = isfinite(r.mean_wind_speed) ? string(round(r.mean_wind_speed; digits=3)) : "NA"
        ri = isfinite(r.mean_ri) ? string(round(r.mean_ri; digits=3)) : "NA"
        write(io, "| $(r.campaign) | $(r.observations) | $(r.height_levels) | $(wind) | $(ri) |\n")
    end
    write(io, "\n")
    write(io, "Overview figure: $(basename(overview_fig))\n\n")
    write(io, "## Output Artifact Manifest\n\n")
    write(io, "| Campaign | Status | Rows | Columns | CSV Stats | JSON | NetCDF | Metrics Figure | Heatmap | K_m Ribbon |\n")
    write(io, "|---|---:|---:|---:|---|---|---|---|---|---|\n")
    for r in eachrow(summary)
        write(io, "| $(r.campaign) | $(r.status) | $(r.n_rows) | $(r.n_columns) | $(basename(r.stats_csv)) | $(basename(r.model_json)) | $(basename(r.netcdf)) | $(basename(r.metrics_fig)) | $(basename(r.heatmap_fig)) | $(basename(r.km_uncertainty_fig)) |\n")
    end
end

tex_path = joinpath(TABLE_DIR, "campaign_summary.tex")
open(tex_path, "w") do io
    println(io, "\\begin{tabular}{lcccc}")
    println(io, "\\toprule")
    println(io, "\\textbf{Campaign} & \\textbf{Observations} & \\textbf{Height Levels} & \\textbf{Mean Wind Speed (\$\\mathrm{m\\,s^{-1}}\$)} & \\textbf{Median Richardson No.} \\\\")
    println(io, "\\midrule")
    for r in eachrow(overview)
        wind = isfinite(r.mean_wind_speed) ? string(round(r.mean_wind_speed; digits = 3)) : "NA"
        ri = isfinite(r.mean_ri) ? string(round(r.mean_ri; digits = 3)) : "NA"
        println(io, "\\texttt{$(r.campaign)} & $(r.observations) & $(r.height_levels) & $(wind) & $(ri) \\\\")
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")
end

println("Campaign export run complete.")
println("Output root: ", OUTPUT_ROOT)
println("Summary CSV: ", summary_csv)
println("Markdown table: ", md_path)
println("LaTeX table: ", tex_path)