"""
    available_ramps()

Return `(:sine_squared, :linear, :smoothstep)`, the built-in normalized ramp
selectors. These named ramps do not accept kwargs. [`ramped_flattop_pulse`](@ref)
also accepts existing pulses and keyword-based pulse constructors.
"""
available_ramps() = (:sine_squared, :linear, :smoothstep)

function _named_ramp(ramp::Symbol, x)
    ramp === :sine_squared && return sinpi(x / 2)^2
    ramp === :linear && return x
    ramp === :smoothstep && return x^2 * (3 - 2x)
    throw(ArgumentError("Unknown ramp $ramp; choose one of $(available_ramps())"))
end

function _flattop_value(p, t, ramp_value)
    (; amplitude, ramp_time, offset, start, stop, duration) = p
    _require_finite("amplitude", amplitude)
    _require_finite("offset", offset)
    _require_finite_real("ramp_time", ramp_time)
    start, stop = _pulse_window(start, stop, duration)
    0 <= 2ramp_time <= stop - start ||
        throw(ArgumentError("Ramps must fit inside the pulse window"))
    start <= t <= stop || return offset

    envelope = if ramp_time == 0
        one(t)
    elseif t < start + ramp_time
        ramp_value((t - start) / ramp_time, true)
    elseif t > stop - ramp_time
        ramp_value((stop - t) / ramp_time, false)
    else
        one(t)
    end
    return offset + amplitude * envelope
end

function _ramped_flattop_pulse(p, t)
    p.ramp in available_ramps() ||
        throw(ArgumentError("Unknown ramp $(p.ramp); choose one of $(available_ramps())"))
    return _flattop_value(p, t, (x, _) -> _named_ramp(p.ramp, x))
end

_register_internal_pulse!(:ramped_flattop, _ramped_flattop_pulse)

function _resolve_ramp(source, kwargs::NamedTuple)
    if source isa Symbol
        source in available_ramps() ||
            throw(ArgumentError("Unknown ramp $source; choose one of $(available_ramps())"))
        isempty(kwargs) || throw(ArgumentError("Named ramps do not accept ramp kwargs"))
        return source
    elseif source isa InternalPulseFunction
        return InternalPulseFunction(source.name, merge(source.parameters, kwargs))
    elseif source isa GenericPulseFunction
        return GenericPulseFunction(source.callable, merge(source.parameters, kwargs))
    elseif source isa AbstractPulse
        isempty(kwargs) || throw(ArgumentError("This ramp pulse does not expose named parameters"))
        return source
    elseif source isa Function
        pulse = source(; merge((; amplitude = 1), kwargs)...)
        pulse isa AbstractPulse || throw(ArgumentError("Ramp factory must return an AbstractPulse"))
        return pulse
    end
    throw(ArgumentError("Ramp must be a pulse, pulse factory, or named ramp"))
end

function _ramp_side_value(source::AbstractPulse, x, split, rising)
    _require_finite_real("split", split)
    0 < split < 1 || throw(ArgumentError("split must lie strictly between 0 and 1"))
    endpoint = rising ? 0 : 1
    baseline = pulse_value(source, endpoint, 1)
    peak = pulse_value(source, split, 1)
    peak != baseline ||
        throw(ArgumentError("Ramp pulse must differ from its endpoint at split"))
    sample = rising ? split * x : 1 - (1 - split) * x
    return (pulse_value(source, sample, 1) - baseline) / (peak - baseline)
end

_ramp_side_value(source::Symbol, x, split, rising) = _named_ramp(source, x)

struct PulseRamps{U,D}
    up::U
    down::D
end

function (r::PulseRamps)(p, t)
    return _flattop_value(p, t, (x, rising) -> rising ?
        _ramp_side_value(r.up, x, p.split_up, true) :
        _ramp_side_value(r.down, x, p.split_down, false))
end

"""
    ramped_flattop_pulse(amplitude, ramp_time; ramp=:sine_squared,
                         ramp_kwargs=(;), ramp_up=ramp, ramp_down=ramp,
                         ramp_up_kwargs=ramp_kwargs, ramp_down_kwargs=ramp_kwargs,
                         split=0.5, split_up=split, split_down=split,
                         offset=0, start=0, stop=nothing)

Create equal-duration ramps around a flat top. Each ramp can be an existing
`AbstractPulse`, a pulse constructor such as `gaussian_pulse` with its own
keyword arguments, or a named shape from `available_ramps()`. Pulse ramps are
sampled over normalized time `0:1`; each `split` marks that pulse's peak. The
rise uses the source pulse before its split; the fall uses the source pulse
after its split. Both are scaled to the requested flat-top amplitude.

Source constructors receive keyword `amplitude=1` plus the supplied NamedTuple
kwargs; existing parameterized pulses receive merged parameter overrides.
Source sigma, center, and frequency use normalized source-time units. Each
pulse split must satisfy `0 < split < 1` and differ in value from the relevant
endpoint. Normalization does not enforce monotonicity or prevent overshoot.

`ramp_time` is the physical duration of each side and must satisfy
`0 ≤ 2ramp_time ≤ stop-start`. Zero ramp time gives an abrupt constant window.
The offset is returned outside the window. The return type is an internal pulse
for matching named ramps with default splits, otherwise a generic pulse.
Resolved source parameters are stored in its callable, not recursively exposed
by [`parameters`](@ref).
"""
function ramped_flattop_pulse(amplitude, ramp_time;
                               ramp = :sine_squared, ramp_kwargs = (;),
                               ramp_up = ramp, ramp_down = ramp,
                               ramp_up_kwargs = ramp_kwargs,
                               ramp_down_kwargs = ramp_kwargs,
                               split = 0.5, split_up = split, split_down = split,
                               offset = 0, start = 0, stop = nothing)
    ramp_up_kwargs isa NamedTuple || throw(ArgumentError("ramp_up_kwargs must be a NamedTuple"))
    ramp_down_kwargs isa NamedTuple || throw(ArgumentError("ramp_down_kwargs must be a NamedTuple"))
    up = _resolve_ramp(ramp_up, ramp_up_kwargs)
    down = _resolve_ramp(ramp_down, ramp_down_kwargs)
    if up isa Symbol && up === down && split_up == 0.5 && split_down == 0.5
        return InternalPulseFunction(:ramped_flattop,
            (; amplitude, ramp_time, offset, start, stop, ramp = up))
    end
    return GenericPulseFunction(PulseRamps(up, down),
        (; amplitude, ramp_time, offset, start, stop, split_up, split_down))
end
