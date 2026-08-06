"""Project normalized observations into a spectral representation (stub)."""
function project_to_gegenbauer(obs::ObservationTable; n_modes::Int = 12, lambda::Float64 = 0.75)
    coeffs = zeros(n_modes)
    return (
        coefficients = coeffs,
        lambda = lambda,
        variables = collect(keys(obs.columns)),
    )
end
