"""Common interface for timed scalar controls, including Pulse and PulseSequence."""
abstract type AbstractPulse end
"""
    Pulse(envelope; duration, amplitude=1, offset=0, delay=0)
    Pulse(envelope, frequency::Number; duration, amplitude=1, offset=0,
          delay=0, phase=0, reference=:pulse)

Envelopes use physical pulse-local time. Active values include the full offset;
outside the closed pulse window the value is offset. A positional frequency mixes
I+iQ as I*cos(θ)-Q*sin(θ). `reference=:sequence` uses time from the immediately
enclosing sequence's onset, including child delay; standalone it uses input time.
EnvelopeSequence duration is inferred and cannot be stretched.
"""
struct Pulse{E,C,A,O,L} <: AbstractPulse
    envelope::E
    carrier::C
    amplitude::A
    offset::O
    delay::L
    function Pulse(e::AbstractEnvelope, c, a, o, d, l)
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
        e = _stretch_envelope(e,d)
        validate_envelope(e, d)
        new{typeof(e),typeof(c),typeof(a),typeof(o),typeof(l)}(e, c, a, o, l)
    end
end
Pulse(e; duration = nothing, amplitude = 1, offset = 0, delay = 0) =
    Pulse(e, nothing, amplitude, offset, duration, delay)
Pulse(e, f::Number; duration = nothing, amplitude = 1, offset = 0, delay = 0, phase = 0, reference = :pulse) =
    Pulse(e, ConstantFrequency(f, phase, reference), amplitude, offset, duration, delay)

Base.getproperty(p::Pulse,k::Symbol)=k===:duration ? _envelope_duration(getfield(p,:envelope)) : getfield(p,k)
Base.propertynames(p::Pulse,private::Bool=false)=(:envelope,:carrier,:amplitude,:offset,:duration,:delay)
