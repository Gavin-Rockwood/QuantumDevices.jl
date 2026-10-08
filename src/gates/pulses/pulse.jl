"""
Timed scalar control interface. Subtypes provide positive `duration`, nonnegative
`delay`, and scalar evaluation as `pulse(t)` with gate-relative time. Built-in
controls use [`Pulse`](@ref); custom waveform shapes normally extend envelopes.
"""
abstract type AbstractPulse end

"""
    Pulse(envelope; duration, delay=0, amplitude=1, offset=0, carrier=nothing)

A timed control, callable as `pulse(t)` with gate-relative time. Its active window
is `[delay, delay + duration]`. Within it, the value is
`offset + amplitude * envelope * carrier`; without a carrier the multiplier is one.
For `IQCarrier`, the carrier mixes the real and imaginary envelope quadratures
into a real signal instead. Outside it the value is `offset`.
Envelopes see time relative to pulse onset.
Duration is positive, delay nonnegative, and both are finite physical times.
Nested envelope/carrier parameters and timing participate in calibration.
"""
struct Pulse{E<:AbstractEnvelope,C,A<:Number,O<:Number,D<:Real,L<:Real} <: AbstractPulse
    envelope::E
    carrier::C
    amplitude::A
    offset::O
    duration::D
    delay::L
    function Pulse(envelope::E, carrier::C, amplitude::A, offset::O, duration::D, delay::L) where {E<:AbstractEnvelope,C,A<:Number,O<:Number,D<:Real,L<:Real}
        carrier isa Union{Nothing,AbstractCarrier} || throw(ArgumentError("carrier must be an AbstractCarrier or nothing"))
        if carrier isa IQCarrier
            _require_finite_real("IQ amplitude", amplitude)
            _require_finite_real("IQ offset", offset)
        end
        _require_finite("amplitude", amplitude)
        _require_finite("offset", offset)
        _require_finite_real("duration", duration) > 0 || throw(ArgumentError("Pulse duration must be positive"))
        _require_finite_real("delay", delay) >= 0 || throw(ArgumentError("Pulse delay must be nonnegative"))
        end_time = _require_finite_real("Pulse end time", delay + duration)
        end_time >= delay || throw(ArgumentError("Pulse end time overflowed"))
        validate_envelope(envelope, duration)
        if carrier !== nothing
            _require_finite("Carrier value", carrier_value(carrier, 0.0, delay))
            _require_finite("Carrier value", carrier_value(carrier, duration, delay + duration))
        end
        pulse = new{E,C,A,O,D,L}(envelope, carrier, amplitude, offset, duration, delay)
        pulse(delay)
        pulse(end_time)
        return pulse
    end
end
Pulse(envelope; duration, delay=0, amplitude=1, offset=0, carrier=nothing) =
    Pulse(envelope, carrier, amplitude, offset, duration, delay)

"""
    pulse_tstops(pulse)
    pulse_tstops(gate::DeviceGate)

Return sorted unique control boundary times for SciML's `tstops` keyword.
Includes pulse onset/end and flattop rise/fall transitions. Gate results contain
only times strictly inside the gate interval. `get_unitary(model, gate)` supplies
these by default; direct `sesolve` calls should pass `tstops=pulse_tstops(gate)`
to avoid skipping short delayed controls. Explicit solver `tstops` override the
default; combine additional times with this list when needed. Custom envelopes
can implement `envelope_tstops` for internal joins or discontinuities.
"""
function pulse_tstops end
