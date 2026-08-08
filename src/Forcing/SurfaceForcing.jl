"""
Observed lower-boundary surface temperature and ground heat flux.

`t_seconds` is elapsed time relative to `t0`, and `Ts` is in Kelvin. Missing
ground heat flux values are represented by `NaN` and are not currently coupled
to the PDE.
"""
struct SurfaceForcing
    t0::DateTime
    t_seconds::Vector{Float64}
    Ts::Vector{Float64}
    G::Vector{Float64}
end

_forcing_datetime(value::DateTime) = value
_forcing_datetime(value) = DateTime(value)

"""
    load_surface_forcing(csv_path, campaign; station=nothing, t0=nothing,
                         source_temperature_unit=:celsius)

Load surface-level observations for one campaign. For station-resolved data,
`station` is required when the campaign contains more than one station.
"""
function load_surface_forcing(
    csv_path::AbstractString,
    campaign::AbstractString;
    station::Union{Nothing,AbstractString}=nothing,
    t0::Union{Nothing,DateTime}=nothing,
    source_temperature_unit::Symbol=:celsius,
)
    isfile(csv_path) || throw(ArgumentError("Surface forcing file not found: $csv_path"))
    data = CSV.read(csv_path, DataFrame)
    required = [:campaign, :timestamp, :G, :Ts]
    missing_columns = setdiff(required, propertynames(data))
    isempty(missing_columns) ||
        throw(ArgumentError("Surface forcing file is missing columns: $(join(missing_columns, ", "))"))

    surface = filter(row -> string(row.campaign) == campaign, data)
    if :level in propertynames(surface)
        filter!(row -> !ismissing(row.level) && string(row.level) == "surface", surface)
    elseif :depth_cm in propertynames(surface)
        filter!(row -> ismissing(row.depth_cm), surface)
    end

    if :station in propertynames(surface)
        available = sort!(unique(string(value) for value in skipmissing(surface.station)))
        selected = station
        if selected === nothing && length(available) > 1
            throw(ArgumentError(
                "Campaign $campaign has multiple stations ($(join(available, ", "))); specify station=...",
            ))
        elseif selected === nothing && length(available) == 1
            selected = only(available)
        end
        selected === nothing || filter!(row -> !ismissing(row.station) && string(row.station) == selected, surface)
    elseif station !== nothing
        throw(ArgumentError("station was provided, but the forcing file has no station column"))
    end

    filter!(row -> !ismissing(row.Ts), surface)
    isempty(surface) && throw(ArgumentError("No surface temperature rows found for campaign=$campaign"))
    surface.timestamp = _forcing_datetime.(surface.timestamp)
    sort!(surface, :timestamp)
    allunique(surface.timestamp) ||
        throw(ArgumentError("Surface forcing timestamps must be unique after campaign/station filtering"))

    origin = t0 === nothing ? first(surface.timestamp) : t0
    elapsed = Dates.value.(surface.timestamp .- origin) ./ 1000.0
    source_temperature_unit in (:celsius, :kelvin) ||
        throw(ArgumentError("source_temperature_unit must be :celsius or :kelvin"))
    temperature_offset = source_temperature_unit === :celsius ? 273.15 : 0.0
    temperatures = Float64.(surface.Ts) .+ temperature_offset
    all(isfinite, temperatures) || throw(ArgumentError("Surface temperatures must be finite"))
    ground_flux = [ismissing(value) ? NaN : Float64(value) for value in surface.G]
    return SurfaceForcing(origin, elapsed, temperatures, ground_flux)
end

"""Linearly interpolate `:Ts` or `:G`, holding endpoint values outside the record."""
function interp_forcing(forcing::SurfaceForcing, field::Symbol, t::Real)
    values = field === :Ts ? forcing.Ts :
             field === :G ? forcing.G :
             throw(ArgumentError("field must be :Ts or :G"))
    times = forcing.t_seconds
    isempty(times) && throw(ArgumentError("Surface forcing cannot be empty"))
    length(times) == length(values) || throw(DimensionMismatch("Forcing time and value lengths differ"))
    length(times) == 1 && return values[1]
    t <= times[1] && return values[1]
    t >= times[end] && return values[end]

    index = clamp(searchsortedlast(times, t), 1, length(times) - 1)
    fraction = (t - times[index]) / (times[index + 1] - times[index])
    return values[index] + fraction * (values[index + 1] - values[index])
end