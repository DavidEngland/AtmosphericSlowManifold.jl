# src/Discovery/WSINDyEngine.jl
struct GegenbauerBasis
    n_spatial::Int
    n_temporal::Int
    lambda::Float64
end

function _discover_feature_key(f::StateVariable)
    return f.name
end

function _discover_feature_key(f::DiagnosticVariable)
    return f.name
end

function _discover_feature_key(f::SpatialDerivative)
    return Symbol("d$(f.order)_$(f.variable)")
end

function _discover_finite_diff_col(col::Vector{Float64}, z::Vector{Float64}, order::Int)
    order >= 1 || return copy(col)
    d = copy(col)
    for _ in 1:order
        out = zeros(length(col))
        out[1] = (d[2] - d[1]) / (z[2] - z[1])
        for i in 2:(length(col) - 1)
            out[i] = (d[i + 1] - d[i - 1]) / (z[i + 1] - z[i - 1])
        end
        out[end] = (d[end] - d[end - 1]) / (z[end] - z[end - 1])
        d = out
    end
    return d
end

function _feature_column_from_obs(obs::ObservationTable, feature::AbstractBasisFeature)
    if feature isa StateVariable || feature isa DiagnosticVariable
        key = _discover_feature_key(feature)
        haskey(obs.columns, key) || throw(KeyError(key))
        return Float64.(obs.columns[key])
    elseif feature isa SpatialDerivative
        haskey(obs.columns, :z) || throw(ArgumentError("ObservationTable must include :z for derivative features."))
        key = feature.variable
        haskey(obs.columns, key) || throw(KeyError(key))
        return _discover_finite_diff_col(Float64.(obs.columns[key]), Float64.(obs.columns[:z]), feature.order)
    end
    throw(ArgumentError("Unsupported basis feature type $(typeof(feature))."))
end

function _solve_bounded_ls(
    G::AbstractMatrix{Float64},
    b::AbstractVector{Float64};
    active::AbstractVector{Int},
    lower_bounds::Union{Nothing, AbstractVector{Float64}} = nothing,
    upper_bounds::Union{Nothing, AbstractVector{Float64}} = nothing,
)
    n_features = size(G, 2)
    model = JuMP.Model(HiGHS.Optimizer)
    JuMP.set_silent(model)

    JuMP.@variable(model, beta[1:n_features])

    inactive = setdiff(collect(1:n_features), collect(active))
    for j in inactive
        JuMP.fix(beta[j], 0.0; force = true)
    end

    if !(lower_bounds === nothing)
        length(lower_bounds) == n_features || throw(ArgumentError("lower_bounds length must match number of features."))
        for j in active
            JuMP.@constraint(model, beta[j] >= lower_bounds[j])
        end
    end

    if !(upper_bounds === nothing)
        length(upper_bounds) == n_features || throw(ArgumentError("upper_bounds length must match number of features."))
        for j in active
            JuMP.@constraint(model, beta[j] <= upper_bounds[j])
        end
    end

    JuMP.@objective(
        model,
        Min,
        0.5 * sum((sum(G[i, j] * beta[j] for j in 1:n_features) - b[i])^2 for i in 1:size(G, 1)),
    )

    JuMP.optimize!(model)
    status = JuMP.termination_status(model)
    status in (JuMP.MOI.OPTIMAL, JuMP.MOI.LOCALLY_SOLVED) ||
        throw(ArgumentError("Constrained least-squares failed with status $(status)."))

    return Float64.(JuMP.value.(beta))
end

"""
    constrained_stlsq(G, b, λ; lower_bounds=nothing, upper_bounds=nothing, max_iter=100, tol=1e-6)

Perform Sequential Thresholded Least Squares with optional box bounds.
Use `lower_bounds` with zeros for non-negativity-constrained terms.
"""
function constrained_stlsq(
    G::AbstractMatrix{Float64},
    b::AbstractVector{Float64},
    λ::Float64;
    lower_bounds::Union{Nothing, AbstractVector{Float64}} = nothing,
    upper_bounds::Union{Nothing, AbstractVector{Float64}} = nothing,
    max_iter::Int = 100,
    tol::Float64 = 1e-6,
)
    size(G, 1) == length(b) || throw(ArgumentError("G and b dimensions are inconsistent."))
    λ >= 0.0 || throw(ArgumentError("λ must be non-negative."))
    max_iter > 0 || throw(ArgumentError("max_iter must be positive."))

    n_features = size(G, 2)
    support = trues(n_features)
    beta = zeros(Float64, n_features)

    for _ in 1:max_iter
        active = findall(identity, support)
        if isempty(active)
            return beta
        end

        beta_new = _solve_bounded_ls(
            Matrix{Float64}(G),
            Vector{Float64}(b);
            active = active,
            lower_bounds = lower_bounds,
            upper_bounds = upper_bounds,
        )

        newly_small = abs.(beta_new) .< λ
        candidate_support = support .& .!newly_small

        if candidate_support == support || norm(beta_new - beta) <= tol
            beta = beta_new
            support = candidate_support
            break
        end

        beta = beta_new
        support = candidate_support
    end

    final_active = findall(identity, support)
    if !isempty(final_active)
        beta = _solve_bounded_ls(
            Matrix{Float64}(G),
            Vector{Float64}(b);
            active = final_active,
            lower_bounds = lower_bounds,
            upper_bounds = upper_bounds,
        )
    end

    beta[abs.(beta) .< λ] .= 0.0
    return beta
end

function _build_feature_evaluation_grid(obs::ObservationTable, library::FeatureLibrary)
    n_features = length(library.features)
    n_features == 0 && return zeros(0, 0)

    cols = [_feature_column_from_obs(obs, f) for f in library.features]
    ns = length(cols[1])
    all(length(c) == ns for c in cols) || throw(ArgumentError("Feature columns in library are inconsistent lengths."))

    grid = zeros(ns, n_features)
    for j in 1:n_features
        grid[:, j] = cols[j]
    end
    return grid
end

"""
Unified Tier-1 discovery pipeline from observation space to `DiscoveredModel`.

Pipeline:
  1) weak-form assembly
  2) physical constraint assembly
  3) sparse optimization
  4) typed OperatorTerm construction
"""
function discover(
    obs::ObservationTable,
    library::FeatureLibrary,
    constraints::Vector{<:AbstractPhysicalConstraint},
    test_family::AbstractTestFunctionFamily,
    optimizer::AbstractSparseOptimizer;
    target_variable::Symbol = :K_m,
    target::Symbol = :u,
    threshold::Float64 = 1e-8,
)
    weak_sys = assemble_weak_system(obs, test_family, library; target = target)

    eval_grid = _build_feature_evaluation_grid(obs, library)
    constraint_matrix = assemble_constraint_matrix(constraints, library.features, eval_grid)

    coeffs = if optimizer isa ConstrainedQP
        solve_sparse_regression(
            weak_sys.G,
            weak_sys.b,
            optimizer;
            A_ineq = constraint_matrix.A_ineq,
            b_ineq = constraint_matrix.b_ineq,
        )
    elseif optimizer isa STRidge
        if !isempty(constraints)
            throw(ArgumentError("STRidge currently does not support inequality constraints; use ConstrainedQP or pass an empty constraint set."))
        end
        solve_sparse_regression(weak_sys.G, weak_sys.b, optimizer)
    else
        throw(ArgumentError("Unsupported optimizer type $(typeof(optimizer))."))
    end

    terms = OperatorTerm{Float64}[]
    for (i, feat) in enumerate(library.features)
        abs(coeffs[i]) <= threshold && continue
        push!(terms, OperatorTerm{Float64}(coeffs[i], BasisOperator[BasisOperator(feat, 1.0)]))
    end

    residual_norm = norm(weak_sys.G * coeffs - weak_sys.b)
    model = DiscoveredModel{Float64}(target_variable, terms, residual_norm, length(terms))
    return model
end

function _candidate_to_basis(expr::Num)
    vars = Symbolics.get_variables(expr)
    if isempty(vars)
        return BasisOperator[]
    end

    if length(vars) == 1
        sym = Symbol(Symbolics.tosymbol(vars[1]))
        return BasisOperator[BasisOperator(StateVariable(sym), 1.0)]
    end

    # Fallback for composite candidates: treat as diagnostic symbolic token.
    label = Symbol(replace(string(expr), r"\s+" => ""))
    return BasisOperator[BasisOperator(DiagnosticVariable(label), 1.0)]
end

"""
Discover a WSINDy closure from observations using the modularized pipeline.

Returns a named tuple with weak matrices, coefficients, discovered IR model,
and a constructed `WSINDyClosure`.
"""
function discover_closure(
    data::ObservationTable,
    candidates::Vector{Num};
    target::Symbol = :u,
    basis::GegenbauerBasis = GegenbauerBasis(),
    lambda::Float64 = 1e-3,
    positivity_constraints::Bool = true,
    diffusivity_indices::Vector{Int} = collect(1:length(candidates)),
    A::Union{Nothing, Matrix{Float64}} = nothing,
    d::Union{Nothing, Vector{Float64}} = nothing,
    threshold::Float64 = 1e-8,
)
    G, b = build_weak_library(data, basis, candidates; target = target)
    coeffs = fit_wsindy_jump(
        G,
        b,
        candidates;
        lambda = lambda,
        positivity_constraints = positivity_constraints,
        diffusivity_indices = diffusivity_indices,
        A = A,
        d = d,
    )

    terms = OperatorTerm{Float64}[]
    for (coef, cand) in zip(coeffs, candidates)
        abs(coef) <= threshold && continue
        push!(terms, OperatorTerm{Float64}(coef, _candidate_to_basis(cand)))
    end

    discovered = DiscoveredModel{Float64}(target, terms, norm(G * coeffs - b), length(terms))
    closure = extract_closure(coeffs, candidates; threshold = threshold)

    return (
        weak_system = WeakFormMatrix(G, b, AbstractBasisFeature[]),
        coefficients = coeffs,
        discovered_model = discovered,
        closure = closure,
    )
end

function GegenbauerBasis(; n_spatial::Int = 8, n_temporal::Int = 1, lambda::Float64 = 0.75)
    n_spatial > 0 || throw(ArgumentError("n_spatial must be positive"))
    n_temporal > 0 || throw(ArgumentError("n_temporal must be positive"))
    return GegenbauerBasis(n_spatial, n_temporal, lambda)
end

"""
Build a weak-form feature matrix G and right-hand side b from observation samples.

This initial engine uses weighted Gegenbauer weak projections in vertical space
and optional polynomial time windows when a `:t` column is present.
"""
function build_weak_library(data::ObservationTable, test_basis::GegenbauerBasis, candidates::Vector{Num}; target::Symbol = :u)
    family = GegenbauerFamily(test_basis.lambda, test_basis.n_spatial)
    wf = assemble_weak_system(data, family, candidates; target = target, n_temporal = test_basis.n_temporal)
    return wf.G, wf.b
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
    A::Union{Nothing, Matrix{Float64}} = nothing,
    d::Union{Nothing, Vector{Float64}} = nothing,
)
    size(G, 1) == length(b) || throw(ArgumentError("G and b dimensions are inconsistent."))
    size(G, 2) == length(candidates) || throw(ArgumentError("G columns must match candidate count."))

    p = length(candidates)
    A_ineq = Matrix{Float64}(undef, 0, p)
    b_ineq = Float64[]

    if positivity_constraints
        A_pos = zeros(length(diffusivity_indices), p)
        for (r, j) in enumerate(diffusivity_indices)
            A_pos[r, j] = 1.0
        end
        A_ineq = vcat(A_ineq, A_pos)
        b_ineq = vcat(b_ineq, zeros(length(diffusivity_indices)))
    end

    if !(A === nothing)
        d === nothing && throw(ArgumentError("d must be provided when A is specified."))
        size(A, 2) == p || throw(ArgumentError("A must have $(p) columns."))
        size(A, 1) == length(d) || throw(ArgumentError("A row count must match length(d)."))
        A_ineq = vcat(A_ineq, A)
        b_ineq = vcat(b_ineq, d)
    end

    solver = ConstrainedQP(lambda = lambda)
    return solve_sparse_regression(G, b, solver; A_ineq = A_ineq, b_ineq = b_ineq)
end

function _candidate_to_state(expr::Num, state::ManifoldState)
    vars = Symbolics.get_variables(expr)
    subs = Dict{Num, Any}()

    for v in vars
        s = lowercase(String(Symbolics.tosymbol(v)))
        if s == "z"
            subs[v] = state.z
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