carrier_value(c::SineCarrier, local_time, gate_time) =
    sinpi(2c.frequency * (c.reference === :pulse ? local_time : gate_time) + c.phase / pi)
carrier_value(c::Carrier, local_time, gate_time) =
    c.callable(c.parameters, c.reference === :pulse ? local_time : gate_time)

function carrier_value(c::IQCarrier, local_time, gate_time)
    clock = c.reference === :pulse ? local_time : gate_time
    angle = 2c.frequency*clock+c.phase/pi
    return complex(sinpi(angle), cospi(angle))
end
_modulate_envelope(::Nothing, shape, local_time, gate_time) = shape
_modulate_envelope(c::AbstractCarrier, shape, local_time, gate_time) =
    shape*_require_finite("Carrier value", carrier_value(c, local_time, gate_time))
function _modulate_envelope(c::IQCarrier, shape, local_time, gate_time)
    quadratures = _require_finite("Carrier value", carrier_value(c, local_time, gate_time))
    return real(shape)*real(quadratures)+imag(shape)*imag(quadratures)
end

function (pulse::Pulse)(t::Real)
    _require_finite_real("time", t)
    pulse.delay <= t <= pulse.delay + pulse.duration || return pulse.offset
    local_time = t - pulse.delay
    # Use the exact stored duration at the end even after floating-point addition.
    t == pulse.delay + pulse.duration && (local_time = pulse.duration)
    shape = pulse.envelope(local_time, pulse.duration)
    signal = _modulate_envelope(pulse.carrier, shape, local_time, t)
    return _require_finite("Pulse value", pulse.offset + pulse.amplitude * signal)
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
    stops = envelope_tstops(pulse.envelope, pulse.duration)
    all(t -> t isa Real && isfinite(t), stops) ||
        throw(ArgumentError("Envelope stops must be finite real times"))
    return sort!(unique(vcat([pulse.delay, pulse.delay+pulse.duration],
        [pulse.delay+t for t in stops if 0 < t < pulse.duration])))
end

Base.show(io::IO, pulse::Pulse) = _show_pulse_summary(io, pulse)

function Base.show(io::IO, ::MIME"text/plain", pulse::Pulse)
    get(io, :compact, false) && return show(io, pulse)
    println(io, "Pulse")
    _show_pulse_details(io, pulse)
end
