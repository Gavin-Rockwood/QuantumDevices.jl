function _gaussian_pulse(p, t)
    (; amplitude, sigma, center, offset, start, stop, duration) = p
    _require_finite("amplitude", amplitude)
    _require_finite("offset", offset)
    _require_finite_real("sigma", sigma)
    sigma > 0 || throw(ArgumentError("Gaussian sigma must be positive"))
    start, stop = _pulse_window(start, stop, duration)
    center = center === nothing ? (start + stop) / 2 : center
    _require_finite_real("center", center)
    start <= center <= stop || throw(ArgumentError("Gaussian center must lie inside its pulse window"))
    start <= t <= stop || return offset
    return offset + amplitude * exp(-(t - center)^2 / (2sigma^2))
end

_register_internal_pulse!(:gaussian, _gaussian_pulse)

"""
    gaussian_pulse(amplitude, sigma; center=nothing, offset=0, start=0, stop=nothing)
    gaussian_pulse(; amplitude=1, sigma, kwargs...)

Return `offset + amplitude * exp(-(t-center)^2/(2sigma^2))` in the inclusive
physical window and `offset` outside it. Positive finite `sigma` uses the same
time units as the gate. Default center is the window midpoint and must lie within
it. The Gaussian is truncated without baseline subtraction, so its edges generally
have a jump. A flattop pulse normalizes source Gaussian endpoints separately.
"""
gaussian_pulse(amplitude, sigma; center = nothing, offset = 0, start = 0, stop = nothing) =
    InternalPulseFunction(:gaussian, (; amplitude, sigma, center, offset, start, stop))

gaussian_pulse(; amplitude = 1, sigma, kwargs...) =
    gaussian_pulse(amplitude, sigma; kwargs...)
