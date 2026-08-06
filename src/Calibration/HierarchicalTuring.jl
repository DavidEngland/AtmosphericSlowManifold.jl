struct CalibrationConfig
    global_scale::Float64
    site_scale::Float64
end

function calibrate_site!(args...; kwargs...)
    throw(ArgumentError("Hierarchical Bayesian calibration is not implemented in MVP scaffold."))
end
