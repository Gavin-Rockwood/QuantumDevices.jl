"""
    AbstractCarrier

Scalar carrier interface. Implement `carrier_value(carrier, local_time, gate_time)`
for a subtype. Frequencies are cycles per unit time and phases are radians.
"""
abstract type AbstractCarrier end

"""Evaluate a carrier with both pulse-local and gate-relative clocks."""
function carrier_value end

function _check_reference(reference)
    reference in (:pulse, :gate) || throw(ArgumentError("Carrier reference must be :pulse or :gate"))
    return reference
end

"""
    SineCarrier(frequency; phase=0, reference=:pulse)

`sin(2π * frequency * clock + phase)`. The clock defaults to time since pulse
onset; `reference=:gate` uses time since gate start. Negative frequencies are valid.
"""
struct SineCarrier{F<:Real,P<:Real} <: AbstractCarrier
    frequency::F
    phase::P
    reference::Symbol
    function SineCarrier(frequency::F, phase::P, reference::Symbol) where {F<:Real,P<:Real}
        _require_finite_real("frequency", frequency)
        _require_finite_real("phase", phase)
        _check_reference(reference)
        new{F,P}(frequency, phase, reference)
    end
end
SineCarrier(frequency; phase=0, reference=:pulse) = SineCarrier(frequency, phase, reference)

"""
    Carrier(f, parameters=(;); reference=:pulse)

Custom scalar carrier. `f(p, t)` receives explicit parameters and the selected
clock (`:pulse` or `:gate`). Callable persistence has the same requirements as
[`Envelope`](@ref). Pulse amplitude and offset stay outside the carrier.
"""
struct Carrier{F,P<:NamedTuple} <: AbstractCarrier
    callable::F
    parameters::P
    reference::Symbol
    function Carrier(f::F, p::P, reference::Symbol) where {F,P<:NamedTuple}
        _check_reference(reference)
        applicable(f, p, 0.0) || throw(ArgumentError("Carrier must accept f(p, t)"))
        new{F,P}(f, p, reference)
    end
end
Carrier(f, parameters=(;); reference=:pulse) = Carrier(f, parameters, reference)
