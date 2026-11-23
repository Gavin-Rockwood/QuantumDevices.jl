"""Common interface for timed scalar controls, including Pulse and PulseSequence."""
abstract type AbstractPulse end
"""Immutable physical definition compiled by `Pulse(recipe)`."""
abstract type AbstractPulseRecipe end
"""
    PulseRecipe(envelope; duration, amplitude=1, offset=0, delay=0, carrier=nothing)
    Pulse(envelope, frequency::Number; duration, amplitude=1, offset=0,
          delay=0, phase=0, reference=:pulse)

Envelopes use physical pulse-local time. Active values include the full offset;
outside the closed pulse window the value is offset. A positional frequency mixes
I+iQ as I*cos(θ)-Q*sin(θ). `reference=:sequence` uses time from the immediately
enclosing sequence's onset, including child delay; standalone it uses input time.
EnvelopeSequence duration is inferred and cannot be stretched.
"""
struct PulseRecipe{E,C,A,O,L} <: AbstractPulseRecipe
    envelope::E
    carrier::C
    amplitude::A
    offset::O
    delay::L
    function PulseRecipe(e::AbstractEnvelope, c, a, o, d, l)
        d === nothing && (d=_envelope_duration(e))
        d === nothing && throw(ArgumentError("duration is required"))
        _require_positive_width("duration", d)
        _require_finite_real("delay", l)>=0 || throw(ArgumentError("delay must be nonnegative"))
        _checked_time_add(d,l)
        _require_finite("amplitude", a)
        _validate_pulse_offset(o)
        c isa Union{Nothing,ConstantFrequency} || throw(ArgumentError("Invalid carrier"))
        if c !== nothing
            _require_finite_real("amplitude", a)
            o isa Num || _require_finite_real("offset", o)
        end
        e = _snapshot_pulse_value(_stretch_envelope(e,d))
        validate_envelope(e, d)
        new{typeof(e),typeof(c),typeof(a),typeof(o),typeof(l)}(e, c, a, o, l)
    end
end
PulseRecipe(e; duration = nothing, amplitude = 1, offset = 0, delay = 0, carrier=nothing) =
    PulseRecipe(e, carrier, amplitude, offset, duration, delay)
PulseRecipe(e, f::Number; duration = nothing, amplitude = 1, offset = 0, delay = 0, phase = 0, reference = :pulse) =
    PulseRecipe(e, ConstantFrequency(f, phase, reference), amplitude, offset, duration, delay)

Base.getproperty(p::PulseRecipe,k::Symbol)=k===:duration ? _envelope_duration(getfield(p,:envelope)) : getfield(p,k)
Base.propertynames(p::PulseRecipe,private::Bool=false)=(:envelope,:carrier,:amplitude,:offset,:duration,:delay)


_snapshot_pulse_value(x) = freeze_recipe_value(x)
_snapshot_pulse_value(e::Envelope) = Envelope(Val(:raw),_snapshot_pulse_value(e.shape),e.duration,e.rise,e.fall)
_snapshot_pulse_value(s::CustomShape) = CustomShape(s.callable,freeze_recipe_value(s.parameters))
_snapshot_pulse_value(s::Region) = Region(_snapshot_pulse_value(s.shape),s.interval)
_snapshot_pulse_value(s::RampedFlattop) = RampedFlattop(Val(:raw),s.rise_time,s.fall_time,_snapshot_pulse_value(s.ramp_up),_snapshot_pulse_value(s.ramp_down))
_snapshot_pulse_value(s::NamedTuple) = map(_snapshot_pulse_value,s)
_snapshot_pulse_value(s::Tuple) = map(_snapshot_pulse_value,s)
_snapshot_pulse_value(e::EnvelopeSequence) = EnvelopeSequence((Segment(_snapshot_pulse_value(s.shape);duration=s.duration,scale=s.scale) for s in e.segments)...)

_snapshot_pulse_value(e::LegacyEnvelope) = e
