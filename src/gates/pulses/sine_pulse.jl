function _sine_pulse(p, t)
    (; amplitude, frequency, phase, offset, start, stop, duration) = p
    _require_finite("amplitude", amplitude)
    _require_finite("offset", offset)
    _require_finite_real("frequency", frequency)
    _require_finite_real("phase", phase)
    start, stop = _pulse_window(start, stop, duration)
    start <= t <= stop || return offset
    return offset + amplitude * sin(2pi * frequency * (t - start) + phase)
end

_register_internal_pulse!(:sine, _sine_pulse)

"""
    sine_pulse(amplitude, frequency; phase=0, offset=0, start=0, stop=nothing)
    sine_pulse(; amplitude=1, frequency, kwargs...)

Return `offset + amplitude * sin(2π * frequency * (t-start) + phase)` within
the inclusive physical window, and offset outside. Frequency is in cycles per
unit time; phase is in radians. This frequency convention differs from the
angular-frequency Hamiltonian convention. All scalar inputs must be finite.
"""
sine_pulse(amplitude, frequency; phase = 0, offset = 0, start = 0, stop = nothing) =
    InternalPulseFunction(:sine, (; amplitude, frequency, phase, offset, start, stop))

sine_pulse(; amplitude = 1, frequency, kwargs...) =
    sine_pulse(amplitude, frequency; kwargs...)
