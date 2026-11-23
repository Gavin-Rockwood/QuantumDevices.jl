"""
    PulseSequence(pulses...; offset=0, delay=0)

Sequential controls. First child delay is measured from sequence onset; each
subsequent delay is a gap after the previous child. Gaps return sequence offset;
active children supply their full value. Outer delay places the whole sequence.
Internal boundaries belong to the following child, except the final endpoint.
"""
struct PulseSequenceRecipe{P,O,D,L} <: AbstractPulseRecipe
    pulses::P
    offset::O
    duration::D
    delay::L
    function PulseSequenceRecipe(children...; offset = 0, delay = 0)
        pulses=map(p -> p isa Pulse ? p.recipe : p,children)
        all(p -> p isa AbstractPulseRecipe,pulses) || throw(ArgumentError("sequence children must be pulse recipes or pulses"))
        isempty(pulses) && throw(ArgumentError("PulseSequence needs children"))
        _validate_pulse_offset(offset)
        _require_finite_real("delay", delay)>=0 || throw(ArgumentError("Negative delay"))
        d=foldl(_checked_time_add,(_checked_time_add(p.delay,p.duration) for p in pulses))
        _checked_time_add(d,delay)
        new{typeof(pulses),typeof(offset),typeof(d),typeof(delay)}(pulses, offset, d, delay)
    end
end

"""Compile sequential children as a `Pulse` with a `PulseSequenceRecipe`."""
PulseSequence(pulses...;kwargs...) = Pulse(PulseSequenceRecipe(pulses...;kwargs...))
