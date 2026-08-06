struct WeakFormMatrix
    G::Matrix{Float64}
    b::Vector{Float64}
    feature_map::Vector{AbstractBasisFeature}
end

function _normalize_axis(v::Vector{Float64})
    vmin, vmax = extrema(v)
    vmax > vmin || throw(ArgumentError("Axis must be non-degenerate for weak projection."))
    return @. 2.0 * (v - vmin) / (vmax - vmin) - 1.0
end

function _trapz(x::Vector{Float64}, y::Vector{Float64})
    length(x) == length(y) || throw(ArgumentError("x and y length mismatch in trapz."))
    acc = 0.0
    for i in 1:(length(x) - 1)
        acc += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i])
    end
    return acc
end

function _wf_feature_key(f::StateVariable)
    return f.name
end

function _wf_feature_key(f::DiagnosticVariable)
    return f.name
end

function _wf_feature_key(f::SpatialDerivative)
    return Symbol("d$(f.order)_$(f.variable)")
end

function _finite_diff_col(col::Vector{Float64}, z::Vector{Float64}, order::Int)
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

function _feature_values(feature::AbstractBasisFeature, obs::ObservationTable)
    haskey(obs.columns, :z) || throw(ArgumentError("ObservationTable must have :z for weak-form assembly."))
    z = obs.columns[:z]

    if feature isa StateVariable || feature isa DiagnosticVariable
        key = _wf_feature_key(feature)
        haskey(obs.columns, key) || throw(KeyError(key))
        return Float64.(obs.columns[key])
    elseif feature isa SpatialDerivative
        key = feature.variable
        haskey(obs.columns, key) || throw(KeyError(key))
        col = Float64.(obs.columns[key])
        return _finite_diff_col(col, z, feature.order)
    else
        throw(ArgumentError("Unsupported feature type $(typeof(feature))"))
    end
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

function _spatial_modes(test_family::AbstractTestFunctionFamily)
    if test_family isa GegenbauerFamily
        return test_family.max_mode
    elseif test_family isa BSplineFamily
        return test_family.num_knots
    end
    throw(ArgumentError("Unsupported test family type $(typeof(test_family))."))
end

function assemble_weak_system(
    obs::ObservationTable,
    test_family::AbstractTestFunctionFamily,
    library::FeatureLibrary;
    target::Symbol = :u,
    n_temporal::Int = 1,
)
    haskey(obs.columns, :z) || throw(ArgumentError("ObservationTable must contain :z for weak system assembly."))
    haskey(obs.columns, target) || throw(ArgumentError("ObservationTable is missing target column $(target)."))

    z = Float64.(obs.columns[:z])
    y = Float64.(obs.columns[target])
    n = length(z)
    n >= 3 || throw(ArgumentError("Need at least 3 vertical levels for weak integration."))

    tvals = haskey(obs.columns, :t) ? Float64.(obs.columns[:t]) : collect(1.0:n)
    z0, H = extrema(z)
    t0, T = extrema(tvals)
    τ = _normalize_axis(tvals)

    n_spatial = _spatial_modes(test_family)
    nrows = n_spatial * n_temporal
    p = length(library.features)

    G = zeros(nrows, p)
    b = zeros(nrows)

    feature_vals = [_feature_values(f, obs) for f in library.features]

    row = 0
    for is in 1:n_spatial
        psi = [evaluate_test_function(test_family, is - 1, zi, z0, H) for zi in z]
        for jt in 1:n_temporal
            row += 1
            omega = @. τ^(jt - 1)
            phi = psi .* omega

            b[row] = _trapz(z, phi .* y)
            for cidx in 1:p
                G[row, cidx] = _trapz(z, phi .* feature_vals[cidx])
            end
        end
    end

    return WeakFormMatrix(G, b, library.features)
end

function assemble_weak_system(
    data::ObservationTable,
    test_family::GegenbauerFamily,
    candidates::Vector{Num};
    target::Symbol = :u,
    n_temporal::Int = 1,
)
    haskey(data.columns, :z) || throw(ArgumentError("ObservationTable must contain :z for weak library construction."))
    haskey(data.columns, target) || throw(ArgumentError("ObservationTable is missing target column $(target)."))
    isempty(candidates) && throw(ArgumentError("candidates cannot be empty."))

    z = Float64.(data.columns[:z])
    y = Float64.(data.columns[target])
    n = length(z)
    n >= 3 || throw(ArgumentError("Need at least 3 vertical levels for weak integration."))

    z0, H = extrema(z)
    tnorm = haskey(data.columns, :t) ? _normalize_axis(Float64.(data.columns[:t])) : ones(n)

    nrows = test_family.max_mode * n_temporal
    p = length(candidates)
    G = zeros(nrows, p)
    b = zeros(nrows)

    row = 0
    for is in 1:test_family.max_mode
        psi = [evaluate_test_function(test_family, is - 1, zi, z0, H) for zi in z]
        for jt in 1:n_temporal
            row += 1
            omega = @. tnorm^(jt - 1)
            phi = psi .* omega

            b[row] = _trapz(z, phi .* y)

            for cidx in 1:p
                vals = [Float64(_candidate_value(candidates[cidx], i, z, data)) for i in 1:n]
                G[row, cidx] = _trapz(z, phi .* vals)
            end
        end
    end

    return WeakFormMatrix(G, b, AbstractBasisFeature[])
end
