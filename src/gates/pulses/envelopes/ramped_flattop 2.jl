const _RAMP_NAMES = (:sine_squared, :linear, :smoothstep)
function _validate_ramp(source, split, rising)
    _require_finite_real("split", split)
    0 < split < 1 || throw(ArgumentError("Ramp split must lie strictly between 0 and 1"))
    source isa Symbol && source in _RAMP_NAMES && return nothing
    source isa AbstractEnvelope || throw(ArgumentError("Ramp must be an envelope or one of $_RAMP_NAMES"))
    validate_envelope(source, 1.0)
    _require_finite("Ramp value", _ramp_value(source, 0.5, split, rising))
    return nothing
end

"""
    RampedFlattop(ramp_time; ramp=:sine_squared, ramp_up=ramp, ramp_down=ramp,
                 rise_time=ramp_time, fall_time=ramp_time,
                 split=0.5, split_up=split, split_down=split)

Unit-height flattop with independent physical rise/fall times. Ramps default to
sin²; `:linear`, `:smoothstep`, and constructed envelope objects are supported.
Reused envelopes are sampled on normalized source time `0:1`; each side subtracts
its endpoint and normalizes its split value to one. Source widths therefore use
normalized units. This does not force monotonicity or prevent overshoot.
All source parameters remain discoverable, including nested ramps.
"""
struct RampedFlattop{U,D,R<:Real,F<:Real,SU<:Real,SD<:Real} <: AbstractEnvelope
    rise_time::R
    fall_time::F
    ramp_up::U
    ramp_down::D
    split_up::SU
    split_down::SD
    function RampedFlattop(rise::R, fall::F, up::U, down::D, su::SU, sd::SD) where {R<:Real,F<:Real,U,D,SU<:Real,SD<:Real}
        _require_finite_real("rise_time", rise) >= 0 || throw(ArgumentError("rise_time must be nonnegative"))
        _require_finite_real("fall_time", fall) >= 0 || throw(ArgumentError("fall_time must be nonnegative"))
        _validate_ramp(up, su, true)
        _validate_ramp(down, sd, false)
        new{U,D,R,F,SU,SD}(rise, fall, up, down, su, sd)
    end
end
function RampedFlattop(ramp_time; ramp=:sine_squared, ramp_up=ramp, ramp_down=ramp,
                      rise_time=ramp_time, fall_time=ramp_time,
                      split=0.5, split_up=split, split_down=split)
    return RampedFlattop(rise_time, fall_time, ramp_up, ramp_down, split_up, split_down)
end

function _ramp_value(source::Symbol, x, split, rising)
    x <= 0 && return zero(x)
    x >= 1 && return one(x)
    source === :sine_squared && return sinpi(x / 2)^2
    source === :linear && return x
    source === :smoothstep && return x^2 * (3 - 2x)
    throw(ArgumentError("Unknown ramp $source"))
end
function _ramp_value(source::AbstractEnvelope, x, split, rising)
    baseline = source(rising ? 0.0 : 1.0, 1.0)
    peak = source(split, 1.0)
    peak != baseline || throw(ArgumentError("Ramp split must differ from its endpoint"))
    x <= 0 && return zero(peak)
    x >= 1 && return one(peak)
    sample = rising ? split * x : 1 - (1 - split) * x
    return (source(sample, 1.0) - baseline) / (peak - baseline)
end
