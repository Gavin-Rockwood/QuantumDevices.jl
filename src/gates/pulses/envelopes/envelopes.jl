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
    value===nothing && throw(ArgumentError("$name is required: construct the shape with duration=T before passing it to Pulse or Segment"))
    _require_finite_real(name, value) > 0 || throw(ArgumentError("$name must be positive"))
    return value
end

function _validated_envelope_tstops(shape::AbstractEnvelope, duration)
    stops = collect(envelope_tstops(shape, duration))
    all(t -> t isa Real && isfinite(t) && 0 <= t <= duration, stops) ||
        throw(ArgumentError("Envelope stops must be finite real times inside its domain"))
    return stops
end

# Checked arithmetic keeps integer timing overflow from silently wrapping.
function _checked_time_add(a,b)
    value=if a isa Base.BitInteger && b isa Base.BitInteger
        x,y=promote(a,b)
        Base.Checked.checked_add(x,y)
    else
        a+b
    end
    return _require_finite_real("time sum",value)
end
function _checked_time_mul(a,b)
    value=if a isa Base.BitInteger && b isa Base.BitInteger
        x,y=promote(a,b)
        Base.Checked.checked_mul(x,y)
    else
        a*b
    end
    return _require_finite_real("time product",value)
end
include("envelope.jl")
include("constant.jl")
include("gaussian.jl")
include("sine_squared.jl")
include("gaussian_zero.jl")
include("sech.jl")
include("cosine.jl")
include("blackman.jl")
include("bump.jl")
include("drag.jl")
include("slepian.jl")
include("gaussian_square.jl")
include("erf_square.jl")
include("InterpolatedEnvelope.jl")
include("EnvelopeSequence.jl")
include("ramped_flattop.jl")
include("overloads.jl")
include("regions.jl")
include("legacy.jl")
