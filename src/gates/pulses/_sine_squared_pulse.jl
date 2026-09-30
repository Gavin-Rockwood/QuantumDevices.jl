function _sine_squared_pulse(p, t)
    (; amplitude, ramp_time, offset, start, stop, duration) = p
    _require_finite("amplitude", amplitude)
    _require_finite("offset", offset)
    _require_finite_real("ramp_time", ramp_time)
    start, stop = _pulse_window(start, stop, duration)
    0 <= 2ramp_time <= stop - start ||
        throw(ArgumentError("Sine-squared ramps must fit inside the pulse window"))
    start <= t <= stop || return offset

    envelope = if ramp_time == 0
        one(t)
    elseif t < start + ramp_time
        sinpi((t - start) / (2ramp_time))^2
    elseif t > stop - ramp_time
        sinpi((stop - t) / (2ramp_time))^2
    else
        one(t)
    end
    return offset + amplitude * envelope
end

_register_internal_pulse!(:sine_squared, _sine_squared_pulse)

sine_squared_pulse(amplitude, ramp_time; offset = 0, start = 0, stop = nothing) =
    InternalPulseFunction(:sine_squared, (; amplitude, ramp_time, offset, start, stop))
