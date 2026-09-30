function _constant_pulse(p, t)
    (; amplitude, offset, start, stop, duration) = p
    _require_finite("amplitude", amplitude)
    _require_finite("offset", offset)
    start, stop = _pulse_window(start, stop, duration)
    return start <= t <= stop ? offset + amplitude : offset
end

_register_internal_pulse!(:constant, _constant_pulse)

"""
    constant_pulse(value; offset=0, start=0, stop=nothing)

Return an internal pulse equal to `offset + value` inside the inclusive physical
time window and `offset` outside it. `stop=nothing` uses gate duration. The window
must satisfy `0 ≤ start ≤ stop ≤ duration`; value and offset must be finite.
"""
constant_pulse(value; offset = 0, start = 0, stop = nothing) =
    InternalPulseFunction(:constant, (; amplitude = value, offset, start, stop))
