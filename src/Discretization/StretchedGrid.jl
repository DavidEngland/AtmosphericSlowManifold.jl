# src/Discretization/StretchedGrid.jl
"""Generate a stretched vertical grid in [0, H] with near-surface clustering."""
function generate_stretched_grid(N::Int, H::Float64, alpha::Float64)
    s = range(0.0, 1.0; length = N)
    return H .* (tanh.(alpha .* s) ./ tanh(alpha))
end
