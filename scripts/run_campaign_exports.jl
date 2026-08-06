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

function build_model_from_row(df::DataFrame)
    coeff_cols = sort(Symbol.(filter(c -> startswith(String(c), "a_"), names(df))))

    if isempty(coeff_cols)
        available = Set(Symbol.(names(df)))
        coeff_cols = [c for c in Symbol[:phi_obs, :zeta, :most_residual, :profile_curvature] if c in available]
    end

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

for (campaign, rel_path) in CAMPAIGN_SOURCES
    slug = slugify(campaign)
    abs_path = normpath(joinpath(pwd(), rel_path))

    if !isfile(abs_path)
        push!(summary, (campaign, abs_path, "missing", 0, 0, "", "", "", "", ""))
        continue
    end

    df = CSV.read(abs_path, DataFrame)
    rename!(df, Symbol.(names(df)))

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
end

summary_csv = joinpath(TABLE_DIR, "campaign_summary.csv")
export_to_csv(summary_csv, summary)

md_path = joinpath(TABLE_DIR, "campaign_summary.md")
open(md_path, "w") do io
    write(io, "# Campaign Production Summary\n\n")
    write(io, "Generated at UTC: $(string(now(UTC)))\n\n")
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
