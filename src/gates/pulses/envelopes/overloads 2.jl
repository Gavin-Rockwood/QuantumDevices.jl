envelope_value(::Constant, t, duration) = one(t)
function envelope_value(shape::Gaussian, t, duration)
    center = shape.center === nothing ? duration / 2 : shape.center
    return exp(-(t - center)^2 / (2shape.sigma^2))
end
function validate_envelope(shape::Gaussian, duration)
    center = shape.center === nothing ? duration / 2 : shape.center
    0 <= center <= duration || throw(ArgumentError("Gaussian center must lie inside the pulse window"))
    return nothing
end
envelope_value(::SineSquared, t, duration) =
    t == 0 || t == duration ? zero(t) : sinpi(t / duration)^2
envelope_value(shape::Envelope, t, duration) = shape.callable(shape.parameters, t, duration)
function validate_envelope(shape::RampedFlattop, duration)
    shape.rise_time + shape.fall_time <= duration || throw(ArgumentError("Rise and fall must fit inside the pulse duration"))
    return nothing
end
function envelope_value(shape::RampedFlattop, t, duration)
    if shape.rise_time > 0 && t < shape.rise_time
        return _ramp_value(shape.ramp_up, t / shape.rise_time, shape.split_up, true)
    elseif shape.fall_time > 0 && t > duration - shape.fall_time
        return _ramp_value(shape.ramp_down, (duration - t) / shape.fall_time, shape.split_down, false)
    end
    return one(t)
end
(shape::AbstractEnvelope)(t, duration) = _require_finite("Envelope value", envelope_value(shape, t, duration))

envelope_tstops(::AbstractEnvelope, duration) = Float64[]
_ramp_tstops(::Symbol, split, rising) = Float64[]
function _ramp_tstops(source::AbstractEnvelope, split, rising)
    stops = envelope_tstops(source, 1.0)
    all(t -> t isa Real && isfinite(t), stops) ||
        throw(ArgumentError("Envelope stops must be finite real times"))
    return rising ? [t/split for t in stops if 0 < t < split] :
        [(1-t)/(1-split) for t in stops if split < t < 1]
end
function envelope_tstops(shape::RampedFlattop, duration)
    stops = Real[shape.rise_time, duration-shape.fall_time]
    if shape.rise_time > 0
        append!(stops, shape.rise_time .* _ramp_tstops(shape.ramp_up, shape.split_up, true))
    end
    if shape.fall_time > 0
        append!(stops, duration .- shape.fall_time .* _ramp_tstops(shape.ramp_down, shape.split_down, false))
    end
    return sort!(unique(filter(t -> 0 < t < duration, stops)))
end
envelope_tstops(shape::Union{GaussianSquare,ErfSquare}, duration) =
    iszero(shape.ramp_time) ? Float64[] : sort!(unique([shape.ramp_time, duration-shape.ramp_time]))
envelope_tstops(shape::Slepian, duration) =
    [duration*j/(shape.samples-1) for j in 1:shape.samples-2]

# expm1 avoids cancellation when a Gaussian is much wider than the pulse.
function _zero_gaussian(sigma, t, duration)
    (t == 0 || t == duration) && return zero(float(t))
    x = (t-duration/2)/sigma
    h = duration/(2sigma)
    !isfinite(h^2) && return exp(-x^2/2)
    denominator = -expm1(-h^2/2)
    if iszero(denominator)
        # Limit as sigma/duration → infinity.
        u = t/duration
        return 4u*(1-u)
    end
    return exp(-x^2/2) * (-expm1((x^2-h^2)/2)) / denominator
end
function _zero_gaussian_derivative(sigma, t, duration)
    h = duration/(2sigma)
    denominator = -expm1(-h^2/2)
    iszero(denominator) && return 4*(1-2t/duration)/duration
    x = (t-duration/2)/sigma
    !isfinite(x) && return zero(float(t))
    return -x*exp(-x^2/2)/(sigma*denominator)
end
envelope_value(shape::GaussianZero, t, duration) = _zero_gaussian(shape.sigma, t, duration)
envelope_value(shape::DRAG, t, duration) = complex(
    _zero_gaussian(shape.sigma, t, duration),
    shape.beta*_zero_gaussian_derivative(shape.sigma, t, duration))
envelope_value(::Cosine, t, duration) =
    t == 0 || t == duration ? zero(float(t)) : sinpi(t/duration)
envelope_value(::Blackman, t, duration) =
    t == 0 || t == duration ? zero(float(t)) : 0.42-0.5cospi(2t/duration)+0.08cospi(4t/duration)
function envelope_value(shape::Sech, t, duration)
    center = shape.center === nothing ? duration/2 : shape.center
    return inv(cosh((t-center)/shape.width))
end
function validate_envelope(shape::Sech, duration)
    center = shape.center === nothing ? duration/2 : shape.center
    0 <= center <= duration || throw(ArgumentError("Sech center must lie inside the pulse window"))
    return nothing
end
function envelope_value(shape::Bump, t, duration)
    center = shape.center === nothing ? duration/2 : shape.center
    x = (t-center)/(duration/2)
    abs(x) >= 1 && return zero(float(x))
    return exp(shape.k * (x^2/(x^2-1)))
end
function validate_envelope(shape::Bump, duration)
    center = shape.center === nothing ? duration/2 : shape.center
    0 <= center <= duration || throw(ArgumentError("Bump center must lie inside the pulse window"))
    return nothing
end
function validate_envelope(shape::Union{GaussianSquare,ErfSquare}, duration)
    2shape.ramp_time <= duration || throw(ArgumentError("Both ramps must fit inside the pulse duration"))
    return nothing
end
function envelope_value(shape::GaussianSquare, t, duration)
    r = shape.ramp_time
    iszero(r) && return one(float(t))
    r <= t <= duration-r && return one(float(t))
    distance = min(t, duration-t)
    distance >= r && return one(float(t))
    # Left half of a centered, baseline-subtracted Gaussian lobe.
    return _zero_gaussian(shape.sigma, distance, 2r)
end
function envelope_value(shape::ErfSquare, t, duration)
    r = shape.ramp_time
    iszero(r) && return one(float(t))
    r <= t <= duration-r && return one(float(t))
    distance = min(t, duration-t)
    distance <= 0 && return zero(float(t))
    distance >= r && return one(float(t))
    h = r/(2sqrt(2)*shape.sigma)
    denominator = erf(h)
    iszero(denominator) && return distance/r
    return (erf((distance-r/2)/(sqrt(2)*shape.sigma))+denominator)/(2denominator)
end
function envelope_value(shape::Slepian, t, duration)
    t <= 0 && return first(shape.weights)
    t >= duration && return last(shape.weights)
    x = (shape.samples-1)*(t/duration)
    index = min(floor(Int, x)+1, shape.samples-1)
    fraction = x-(index-1)
    return muladd(fraction, shape.weights[index+1]-shape.weights[index], shape.weights[index])
end
parameters(shape::Slepian) = Dict{String,Tuple}(
    "time_bandwidth" => ("time_bandwidth", shape.time_bandwidth))
function _replace_property(shape::Slepian, key::Symbol, value)
    key === :time_bandwidth && return Slepian(value; samples=shape.samples)
    key === :samples && return Slepian(shape.time_bandwidth; samples=value)
    throw(ArgumentError("Slepian.$key is derived or unknown"))
end
function _show_control_shape(io, shape::Slepian)
    print(io, "Slepian(time_bandwidth=", shape.time_bandwidth, ", samples=", shape.samples, ")")
end

Base.show(io::IO, shape::Slepian) = _show_control_shape(io, shape)
Base.show(io::IO, ::MIME"text/plain", shape::Slepian) = _show_control_shape(io, shape)
