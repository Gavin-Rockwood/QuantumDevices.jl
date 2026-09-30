function _sine_squared_pulse(p, t)
    (; amplitude, frequency, phase, offset, start, stop, duration) = p
    _require_finite("amplitude", amplitude)
    _require_finite("offset", offset)
    _require_finite_real("frequency", frequency)
    _require_finite_real("phase", phase)
    start, stop = _pulse_window(start, stop, duration)
    start <= t <= stop || return offset
    return offset + amplitude * sin(2pi * frequency * (t - start) + phase)^2
end

_register_internal_pulse!(:sine_squared, _sine_squared_pulse)

"""
    sine_squared_pulse(amplitude, frequency; phase=0, offset=0, start=0, stop=nothing)
    sine_squared_pulse(; amplitude=1, frequency, kwargs...)

Return `offset + amplitude * sin(2π * frequency * (t-start) + phase)^2` within
the inclusive physical window, and offset outside. Frequency describes the
unsquared sine in cycles per time; its square has twice that repetition rate.
For a single normalized lobe use `frequency=0.5`, `phase=0`. This constructor is
periodic; use [`ramped_flattop_pulse`](@ref) for a flat top.
"""
sine_squared_pulse(amplitude, frequency; phase = 0, offset = 0, start = 0, stop = nothing) =
    InternalPulseFunction(:sine_squared, (; amplitude, frequency, phase, offset, start, stop))

sine_squared_pulse(; amplitude = 1, frequency, kwargs...) =
    sine_squared_pulse(amplitude, frequency; kwargs...)
