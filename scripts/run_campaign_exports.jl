using AtmosphericSlowManifold
using CSV
using DataFrames
using Dates
using Plots
using Statistics

const CAMPAIGN_SOURCES = Dict(
    "CASES-99" => "../SpectralBL-Analytics/data/drafts/trajectories/trajectory_cases_99.csv",
    "FLOSS" => "../SpectralBL-Analytics/data/drafts/trajectories/trajectory_floss.csv",
    "BLLAST" => "../SpectralBL-Analytics/data/drafts/trajectories/trajectory_bllast.csv",
    "SHEBA" => "../SpectralBL-Analytics/data/sheba/processed/sheba_input.csv",
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
        return (n = 0, mean = NaN, std = NaN, min = NaN, max = NaN)
    end
    return (
        n = length(vf),
        mean = mean(vf),
        std = length(vf) > 1 ? std(vf) : 0.0,
        min = minimum(vf),
        max = maximum(vf),
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

function ensure_obukhov_scaling!(df::DataFrame)
    has_l = hasproperty(df, :L_obukhov)
    if has_l
        lvals = get_numeric_column(df, :L_obukhov)
        finite_count(lvals) > 0 && return
    end

    hasproperty(df, :eta_3) || return
    eta3 = get_numeric_column(df, :eta_3)
    zref = choose_reference_height(df)
    L = fill(NaN, length(eta3))
    for i in eachindex(eta3)
        η = eta3[i]
        isfinite(η) || continue
        abs(η) <= 1e-8 && continue
        L[i] = zref / η
    end
    df[!, :L_obukhov] = L
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
        zeta[i] = zref / Li
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
            z_col = :sample_index,
            u_col = :eta_1,
            v_col = :eta_2,
            temp_col = temp_candidate,
            ustar_col = choose_optional_finite_column(df, Symbol[:ustar, :u_star]),
            hs_col = get_existing_column(df, Symbol[:hs, :h, :shf, :sensible_heat_flux]),
        )
    elseif hasproperty(df, :z_lo) && hasproperty(df, :ws_lo) && hasproperty(df, :ws_hi) && hasproperty(df, :T_lo)
        (
            z_col = :z_lo,
            u_col = :ws_lo,
            v_col = :ws_hi,
            temp_col = :T_lo,
            q_col = choose_optional_finite_column(df, Symbol[:q_lo, :q, :q_hi, :specific_humidity]),
            ustar_col = choose_optional_finite_column(df, Symbol[:ustar, :u_star]),
            hs_col = get_existing_column(df, Symbol[:hs, :h, :shf, :sensible_heat_flux]),
            theta_ref_col = get_existing_column(df, Symbol[:T_lo, :theta, :theta_k]),
        )
    else
        (
            z_col = :z,
            u_col = :u,
            v_col = :v,
            temp_col = :theta,
            q_col = choose_optional_finite_column(df, Symbol[:q]),
            ustar_col = choose_optional_finite_column(df, Symbol[:u_star, :ustar]),
            hs_col = get_existing_column(df, Symbol[:hs, :h, :shf, :sensible_heat_flux]),
            theta_ref_col = get_existing_column(df, Symbol[:theta, :theta_k]),
        )
    end

    obs = read_observation_data(
        abs_path;
        kwargs...,
        compute_obukhov = true,
        surface_flux_aliases = true,
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

function collect_stats(df::DataFrame, cols::Vector{Symbol})
    rows = NamedTuple[]
    diag = Dict{Symbol, Any}()
    for c in cols
        v = get_numeric_column(df, c)
        isempty(v) && continue
        st = finite_stats(v)
        push!(rows, (metric = String(c), n = st.n, mean = st.mean, std = st.std, min = st.min, max = st.max))
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
    ri_cols = sort(Symbol.(filter(n -> startswith(String(n), "ri_g_"), names(df))))
    isempty(ri_cols) && return nothing

    z_grid = Float64[]
    for c in ri_cols
        suffix = replace(String(c)[6:end], "_" => ".")
        z = tryparse(Float64, suffix)
        push!(z_grid, z === nothing ? NaN : z)
    end
    if any(!isfinite, z_grid)
        return nothing
    end

    order = sortperm(z_grid)
    z_grid = z_grid[order]
    ri_cols = ri_cols[order]

    t_grid = collect(1.0:1.0:nrow(df))
    mat = zeros(Float64, length(z_grid), nrow(df))
    for (i, c) in enumerate(ri_cols)
        vec = get_numeric_column(df, c)
        length(vec) == nrow(df) || return nothing
        mat[i, :] = vec
    end
    return z_grid, t_grid, mat
end

function ri_height_levels(df::DataFrame)
    rm = ri_matrix(df)
    return rm === nothing ? 0 : length(rm[1])
end

function mean_ri_value(df::DataFrame)
    rm = ri_matrix(df)
    if !(rm === nothing)
        _, _, mat = rm
        return finite_mean(vec(mat))
    elseif hasproperty(df, :zeta)
        return finite_mean(get_numeric_column(df, :zeta))
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
    p = plot(size = (1100, 700), xlabel = "Sample Index", ylabel = "Metric Value", title = "$(campaign): Diagnostic Metrics")
    for c in cols
        y = get_numeric_column(df, c)
        plot!(p, x, y, lw = 2, label = String(c))
    end
    out = joinpath(FIG_DIR, "$(slug)_metrics.png")
    savefig(p, out)
    return out
end

function plot_ri_heatmap(campaign::String, slug::String, z_grid::Vector{Float64}, t_grid::Vector{Float64}, mat::Matrix{Float64})
    p = heatmap(
        t_grid,
        z_grid,
        mat,
        xlabel = "Sample Index",
        ylabel = "Height Proxy z (m)",
        colorbar_title = "Ri_g",
        title = "$(campaign): Ri_g(z, t) Heatmap",
        size = (1100, 700),
    )
    out = joinpath(FIG_DIR, "$(slug)_ri_heatmap.png")
    savefig(p, out)
    return out
end

function plot_campaign_overview(overview::DataFrame)
    campaigns = String.(overview.campaign)
    wind = Float64.(overview.mean_wind_speed)
    ri = Float64.(overview.mean_ri)

    p1 = bar(
        campaigns,
        wind,
        xlabel = "Campaign",
        ylabel = "Mean Wind Speed (m s^-1)",
        title = "Campaign Mean Wind Speed",
        legend = false,
        size = (1100, 450),
        color = :steelblue,
    )

    p2 = bar(
        campaigns,
        ri,
        xlabel = "Campaign",
        ylabel = "Mean Richardson Number",
        title = "Campaign Mean Stability",
        legend = false,
        size = (1100, 450),
        color = :darkorange,
    )

    p = plot(p1, p2; layout = (2, 1), size = (1100, 900))
    out = joinpath(FIG_DIR, "campaign_overview.png")
    savefig(p, out)
    return out
end

summary = DataFrame(
    campaign = String[],
    source_file = String[],
    status = String[],
    n_rows = Int[],
    n_columns = Int[],
    stats_csv = String[],
    model_json = String[],
    netcdf = String[],
    metrics_fig = String[],
    heatmap_fig = String[],
)

overview = DataFrame(
    campaign = String[],
    observations = Int[],
    height_levels = Int[],
    mean_wind_speed = Float64[],
    mean_ri = Float64[],
)

for (campaign, rel_path) in CAMPAIGN_SOURCES
    slug = slugify(campaign)
    abs_path = normpath(joinpath(pwd(), rel_path))

    if !isfile(abs_path)
        push!(summary, (campaign, abs_path, "missing", 0, 0, "", "", "", "", ""))
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

    stats_df, diag_stats = collect_stats(df, Symbol[:most_residual, :profile_curvature, :transversality, :phi_obs, :zeta, :L_obukhov, :ustar])
    stats_csv = joinpath(CSV_DIR, "$(slug)_stats.csv")
    export_to_csv(stats_csv, stats_df)

    model = build_model_from_row(df)
    diagnostics = Dict{Symbol, Any}(
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
    write(io, "| Campaign | Observations | Height Levels | Mean Wind Speed (m s^-1) | Mean Richardson Number |\n")
    write(io, "|---|---:|---:|---:|---:|\n")
    for r in eachrow(overview)
        wind = isfinite(r.mean_wind_speed) ? string(round(r.mean_wind_speed; digits = 3)) : "NA"
        ri = isfinite(r.mean_ri) ? string(round(r.mean_ri; digits = 3)) : "NA"
        write(io, "| $(r.campaign) | $(r.observations) | $(r.height_levels) | $(wind) | $(ri) |\n")
    end
    write(io, "\n")
    write(io, "Overview figure: $(basename(overview_fig))\n\n")
    write(io, "## Output Artifact Manifest\n\n")
    write(io, "| Campaign | Status | Rows | Columns | CSV Stats | JSON | NetCDF | Metrics Figure | Heatmap |\n")
    write(io, "|---|---:|---:|---:|---|---|---|---|---|\n")
    for r in eachrow(summary)
        write(io, "| $(r.campaign) | $(r.status) | $(r.n_rows) | $(r.n_columns) | $(basename(r.stats_csv)) | $(basename(r.model_json)) | $(basename(r.netcdf)) | $(basename(r.metrics_fig)) | $(basename(r.heatmap_fig)) |\n")
    end
end

tex_path = joinpath(TABLE_DIR, "campaign_summary.tex")
open(tex_path, "w") do io
    write(io, "\\begin{table}[htbp]\n")
    write(io, "  \\centering\n")
    write(io, "  \\caption{Atmospheric Boundary Layer Field Campaign Dataset Overview}\n")
    write(io, "  \\label{tab:campaign-overview}\n")
    write(io, "  \\begin{tabular}{lcccc}\n")
    write(io, "    \\toprule\n")
    write(io, "    \\textbf{Campaign} & \\textbf{Observations} & \\textbf{Height Levels} & \\textbf{Mean Wind Speed (\$\\mathrm{m\\,s^{-1}}\$)} & \\textbf{Mean Richardson No.} \\\\ \n")
    write(io, "    \\midrule\n")
    for r in eachrow(overview)
        wind = isfinite(r.mean_wind_speed) ? string(round(r.mean_wind_speed; digits = 3)) : "NA"
        ri = isfinite(r.mean_ri) ? string(round(r.mean_ri; digits = 3)) : "NA"
        write(io, "    \\texttt{$(r.campaign)} & $(r.observations) & $(r.height_levels) & $(wind) & $(ri) \\\\ \n")
    end
    write(io, "    \\bottomrule\n")
    write(io, "  \\end{tabular}\n")
    write(io, "\\end{table}\n\n")

    write(io, "\\begin{figure}[htbp]\n")
    write(io, "  \\centering\n")
    write(io, "  \\includegraphics[width=0.92\\linewidth]{reports/generated/campaign_exports/figures/$(basename(overview_fig))}\n")
    write(io, "  \\caption{Comparative campaign overview for derived mean wind speed and mean stability metrics.}\n")
    write(io, "  \\label{fig:campaign-overview}\n")
    write(io, "\\end{figure}\n\n")

    write(io, "\\begin{table}[htbp]\n")
    write(io, "  \\centering\n")
    write(io, "  \\caption{Campaign Production Output Summary}\n")
    write(io, "  \\label{tab:campaign-production-summary}\n")
    write(io, "  \\begin{tabular}{lccc}\n")
    write(io, "    \\toprule\n")
    write(io, "    \\textbf{Campaign} & \\textbf{Status} & \\textbf{Rows} & \\textbf{Columns} \\\\ \n")
    write(io, "    \\midrule\n")
    for r in eachrow(summary)
        write(io, "    \\texttt{$(r.campaign)} & \\texttt{$(r.status)} & $(r.n_rows) & $(r.n_columns) \\\\ \n")
    end
    write(io, "    \\bottomrule\n")
    write(io, "  \\end{tabular}\n")
    write(io, "\\end{table}\n")
end

println("Campaign export run complete.")
println("Output root: ", OUTPUT_ROOT)
println("Summary CSV: ", summary_csv)
println("Markdown table: ", md_path)
println("LaTeX table: ", tex_path)
