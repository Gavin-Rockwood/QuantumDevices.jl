"""
    AbstractEnvelope

Reusable dimensionless shape evaluated by `envelope_value(shape, t, duration)`.
Time is relative to pulse onset. Implement that method for custom subtypes;
optionally implement `validate_envelope(shape, duration)` and persistence hooks.
"""
abstract type AbstractEnvelope end

"""Evaluate a dimensionless envelope at local time `t` over `duration`."""
function envelope_value end

"""Validate an envelope for a positive pulse duration."""
function validate_envelope(shape::AbstractEnvelope, duration)
    _require_finite("Envelope value", envelope_value(shape, 0.0, duration))
    _require_finite("Envelope value", envelope_value(shape, duration, duration))
    return nothing
end

"""
    envelope_tstops(shape::AbstractEnvelope, duration)

Interior pulse-local times at which adaptive integration should stop. Defaults
to no stops. Custom envelopes can implement this for joins or discontinuities;
`pulse_tstops` adds pulse boundaries and shifts these times by the pulse delay.
"""
function envelope_tstops end

function _require_positive_width(name, value)
    _require_finite_real(name, value) > 0 || throw(ArgumentError("$name must be positive"))
    return value
end

include("constant.jl")
include("gaussian.jl")
include("sine_squared.jl")
include("envelope.jl")
include("ramped_flattop.jl")
include("gaussian_zero.jl")
include("sech.jl")
include("cosine.jl")
include("blackman.jl")
include("bump.jl")
include("drag.jl")
include("slepian.jl")
include("gaussian_square.jl")
include("erf_square.jl")
include("overloads.jl")
