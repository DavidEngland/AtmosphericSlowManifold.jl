using HiGHS

struct GegenbauerBasis
    n_spatial::Int
    n_temporal::Int
    lambda::Float64
end

function GegenbauerBasis(; n_spatial::Int = 8, n_temporal::Int = 1, lambda::Float64 = 0.75)
    n_spatial > 0 || throw(ArgumentError("n_spatial must be positive"))
    n_temporal > 0 || throw(ArgumentError("n_temporal must be positive"))
    return GegenbauerBasis(n_spatial, n_temporal, lambda)
end

function _engine_gegenbauerC(n::Int, lambda::Float64, x::Float64)
    n == 0 && return 1.0
    n == 1 && return 2.0 * lambda * x

    c_nm2 = 1.0
    c_nm1 = 2.0 * lambda * x
    for k in 2:n
        c_n = (2.0 * (k + lambda - 1.0) * x * c_nm1 - (k + 2.0 * lambda - 2.0) * c_nm2) / k
        c_nm2 = c_nm1
        c_nm1 = c_n
    end
    return c_nm1
end

function _normalize_axis(v::Vector{Float64})
    vmin, vmax = extrema(v)
    vmax > vmin || throw(ArgumentError("Axis must be non-degenerate for weak projection."))
    return @. 2.0 * (v - vmin) / (vmax - vmin) - 1.0
end

function _trapz_weighted(x::Vector{Float64}, f::Vector{Float64}, g::Vector{Float64}, lambda::Float64)
    acc = 0.0
    for i in 1:(length(x) - 1)
        x1 = x[i]
        x2 = x[i + 1]
        w1 = max(1e-12, (1.0 - x1^2)^(lambda - 0.5))
        w2 = max(1e-12, (1.0 - x2^2)^(lambda - 0.5))
        y1 = f[i] * g[i] * w1
        y2 = f[i + 1] * g[i + 1] * w2
        acc += 0.5 * (y1 + y2) * (x2 - x1)
    end
    return acc
end

function _candidate_value(expr::Num, i::Int, z::Vector{Float64}, data::ObservationTable)
    vars = Symbolics.get_variables(expr)
    if isempty(vars)
        raw = Symbolics.value(expr)
        if raw isa Number
            return float(raw)
        end
        parsed = tryparse(Float64, string(raw))
        parsed === nothing && throw(ArgumentError("Constant candidate expression could not be evaluated numerically."))
        return parsed
    end

    subs = Dict{Num, Any}()
    for v in vars
        s = lowercase(String(Symbolics.tosymbol(v)))
        if s == "z"
            subs[v] = z[i]
            continue
        end

        # Accept explicit column names when candidates use atmospheric symbols.
        if haskey(data.columns, Symbol(s))
            subs[v] = data.columns[Symbol(s)][i]
            continue
        end

        if s == "u_star" && haskey(data.columns, :u_star)
            subs[v] = data.columns[:u_star][i]
            continue
        end
    end

    val = Symbolics.value(Symbolics.substitute(expr, subs))
    val isa Number || throw(ArgumentError("Candidate expression could not be evaluated numerically at sample index $(i)."))
    return float(val)
end

"""
Build a weak-form feature matrix G and right-hand side b from observation samples.

This initial engine uses weighted Gegenbauer weak projections in vertical space
and optional polynomial time windows when a `:t` column is present.
"""
function build_weak_library(data::ObservationTable, test_basis::GegenbauerBasis, candidates::Vector{Num}; target::Symbol = :u)
    haskey(data.columns, :z) || throw(ArgumentError("ObservationTable must contain :z for weak library construction."))
    haskey(data.columns, target) || throw(ArgumentError("ObservationTable is missing target column $(target)."))
    isempty(candidates) && throw(ArgumentError("candidates cannot be empty."))

    z = data.columns[:z]
    y = data.columns[target]
    n = length(z)
    n >= 3 || throw(ArgumentError("Need at least 3 vertical levels for weak integration."))

    x = _normalize_axis(z)
    tnorm = haskey(data.columns, :t) ? _normalize_axis(data.columns[:t]) : ones(n)

    nrows = test_basis.n_spatial * test_basis.n_temporal
    p = length(candidates)
    G = zeros(nrows, p)
    b = zeros(nrows)

    row = 0
    for is in 1:test_basis.n_spatial
        psi = [_engine_gegenbauerC(is - 1, test_basis.lambda, xv) for xv in x]
        for jt in 1:test_basis.n_temporal
            row += 1
            omega = @. tnorm^(jt - 1)
            phi = psi .* omega

            b[row] = _trapz_weighted(x, phi, y, test_basis.lambda)

            for cidx in 1:p
                vals = [Float64(_candidate_value(candidates[cidx], i, z, data)) for i in 1:n]
                G[row, cidx] = _trapz_weighted(x, phi, vals, test_basis.lambda)
            end
        end
    end

    return G, b
end

"""
Solve a constrained sparse regression problem for WSINDy coefficients.

Objective:
    min 0.5||G*Xi - b||_2^2 + lambda*||Xi||_1
"""
function fit_wsindy_jump(
    G::Matrix{Float64},
    b::Vector{Float64},
    candidates::Vector{Num};
    lambda::Float64 = 1e-3,
    positivity_constraints::Bool = true,
    diffusivity_indices::Vector{Int} = collect(1:length(candidates)),
)
    size(G, 1) == length(b) || throw(ArgumentError("G and b dimensions are inconsistent."))
    size(G, 2) == length(candidates) || throw(ArgumentError("G columns must match candidate count."))

    p = length(candidates)
    model = JuMP.Model(HiGHS.Optimizer)
    JuMP.set_silent(model)

    JuMP.@variable(model, xi[1:p])
    JuMP.@variable(model, t[1:p] >= 0)
    JuMP.@constraint(model, [j in 1:p], t[j] >= xi[j])
    JuMP.@constraint(model, [j in 1:p], t[j] >= -xi[j])

    if positivity_constraints
        for j in diffusivity_indices
            JuMP.@constraint(model, xi[j] >= 0)
        end
    end

    JuMP.@objective(
        model,
        Min,
        0.5 * sum((sum(G[i, j] * xi[j] for j in 1:p) - b[i])^2 for i in 1:size(G, 1)) + lambda * sum(t)
    )

    JuMP.optimize!(model)
    status = JuMP.termination_status(model)
    status in (JuMP.MOI.OPTIMAL, JuMP.MOI.LOCALLY_SOLVED) || throw(ArgumentError("WSINDy optimization failed with status $(status)."))

    return Float64.(JuMP.value.(xi))
end

function _candidate_to_state(expr::Num, state::ManifoldState)
    vars = Symbolics.get_variables(expr)
    subs = Dict{Num, Any}()

    for v in vars
        s = lowercase(String(Symbolics.tosymbol(v)))
        if s == "z"
            subs[v] = state.z0
        elseif s == "r"
            subs[v] = state.r
        elseif s == "chi"
            subs[v] = state.chi
        elseif s == "pi_g"
            subs[v] = state.pi_g
        elseif s == "eta1"
            subs[v] = state.eta1
        elseif s == "eta2"
            subs[v] = state.eta2
        elseif s == "eta3"
            subs[v] = state.eta3
        elseif s == "u_star"
            subs[v] = state.u_star
        end
    end

    return Symbolics.substitute(expr, subs)
end

"""Convert sparse coefficient vector and candidate basis into a WSINDyClosure."""
function extract_closure(coefficients::Vector{Float64}, candidates::Vector{Num}; threshold::Float64 = 1e-8)
    length(coefficients) == length(candidates) || throw(ArgumentError("coefficients and candidates lengths must match."))

    state = ManifoldState()
    km_expr = zero(state.z0)
    for (coef, cand) in zip(coefficients, candidates)
        abs(coef) <= threshold && continue
        km_expr += coef * _candidate_to_state(cand, state)
    end

    kh_expr = km_expr / 0.74
    flux_expr = state.u_star^2
    return WSINDyClosure(km_expr, kh_expr, flux_expr)
end