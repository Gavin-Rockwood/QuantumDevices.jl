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
carrier_value(c::SineCarrier, local_time, gate_time) =
    sinpi(2c.frequency * (c.reference === :pulse ? local_time : gate_time) + c.phase / pi)
carrier_value(c::Carrier, local_time, gate_time) =
    c.callable(c.parameters, c.reference === :pulse ? local_time : gate_time)

function (pulse::Pulse)(t::Real)
    _require_finite_real("time", t)
    pulse.delay <= t <= pulse.delay + pulse.duration || return pulse.offset
    local_time = t - pulse.delay
    # Use the exact stored duration at the end even after floating-point addition.
    t == pulse.delay + pulse.duration && (local_time = pulse.duration)
    shape = pulse.envelope(local_time, pulse.duration)
    carrier = pulse.carrier === nothing ? 1 :
        _require_finite("Carrier value", carrier_value(pulse.carrier, local_time, t))
    return _require_finite("Pulse value", pulse.offset + pulse.amplitude * shape * carrier)
end

function _add_control_parameter!(result, name, path, value)
    if value isa Union{AbstractEnvelope,AbstractCarrier,AbstractPulse}
        for (child_name, (child_path, leaf)) in parameters(value)
            result["$name/$child_name"] = ("$path/$child_path", leaf)
        end
    elseif value isa NamedTuple
        for (key, entry) in pairs(value)
            _add_control_parameter!(result, "$name/$key", "$path/$key", entry)
        end
    else
        result[name] = (path, value)
    end
    return result
end
function _collect_control_parameters!(result, object)
    for key in fieldnames(typeof(object))
        key === :callable && continue
        value = getfield(object, key)
        if key === :parameters && value isa NamedTuple
            for (name, entry) in pairs(value)
                # Keep explicit paths for parameters shadowed by structural fields.
                path = hasproperty(object, name) ? "parameters/$name" : String(name)
                _add_control_parameter!(result, path, path, entry)
            end
        else
            _add_control_parameter!(result, String(key), String(key), value)
        end
    end
    return result
end
parameters(object::Union{AbstractEnvelope,AbstractCarrier,AbstractPulse}) =
    _collect_control_parameters!(Dict{String,Tuple}(), object)

function _validate_timed_control(pulse::AbstractPulse)
    _require_finite_real("Pulse duration", pulse.duration) > 0 || throw(ArgumentError("Pulse duration must be positive"))
    _require_finite_real("Pulse delay", pulse.delay) >= 0 || throw(ArgumentError("Pulse delay must be nonnegative"))
    end_time = _require_finite_real("Pulse end time", pulse.delay + pulse.duration)
    end_time >= pulse.delay || throw(ArgumentError("Pulse end time overflowed"))
    _require_finite("Pulse value", pulse(0.0))
    _require_finite("Pulse value", pulse(pulse.delay))
    _require_finite("Pulse value", pulse(end_time))
    return end_time
end

pulse_tstops(pulse::AbstractPulse) = [pulse.delay, pulse.delay + pulse.duration]
function pulse_tstops(pulse::Pulse)
    shape = pulse.envelope
    if shape isa RampedFlattop
        return sort!(unique([pulse.delay, pulse.delay + shape.rise_time,
            pulse.delay + pulse.duration - shape.fall_time, pulse.delay + pulse.duration]))
    end
    return [pulse.delay, pulse.delay + pulse.duration]
end
