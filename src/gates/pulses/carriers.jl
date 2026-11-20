"""Constant cosine carrier; frequency is cycles/time, phase is radians."""
struct ConstantFrequency{F<:Real,P<:Real}
    frequency::F
    phase::P
    reference::Symbol
    function ConstantFrequency(frequency, phase, reference)
        _require_finite_real("frequency", frequency)
        _require_finite_real("phase", phase)
        reference in (:pulse, :sequence) || throw(ArgumentError("reference must be :pulse or :sequence"))
        new{typeof(frequency),typeof(phase)}(frequency, phase, reference)
    end
end
function _modulate(c::ConstantFrequency, z, local_time, sequence_time)
    t = c.reference === :pulse ? local_time : sequence_time
    angle = 2pi*c.frequency*t+c.phase
    return real(z)*cos(angle)-imag(z)*sin(angle)
end
_modulate(::Nothing, z, local_time, sequence_time) = z
