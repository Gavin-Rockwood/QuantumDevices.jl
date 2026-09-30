"""
    AbstractPulse

Interface for time-dependent scalar gate parameters. Implement
[`pulse_value`](@ref) for a subtype; optionally specialize `validate_pulse`,
[`parameters`](@ref), and persistence hooks. Pulses are callable as
`pulse(t, duration)`. Returned values must be finite scalar numbers.
"""
abstract type AbstractPulse end

"""
    AbstractParameterizedPulse <: AbstractPulse

Pulse interface with a named `parameters` field and [`pulse_function`](@ref).
The shared evaluator merges the owning gate's `duration` into those parameters
and calls `f(p, t)`. Stored parameters must not contain `duration`.
"""
abstract type AbstractParameterizedPulse <: AbstractPulse end

"""
    pulse_function(pulse::AbstractParameterizedPulse)

Resolve the `f(p, t)` callable behind a parameterized pulse. Its first argument
contains the pulse parameters plus the owning gate duration. Internal pulses
resolve a registered name; generic pulses return their stored callable.
"""
function pulse_function end

function _check_pulse_parameters(parameters::NamedTuple)
    haskey(parameters, :duration) && throw(ArgumentError(
        "duration is supplied by the owning gate, not stored in pulse parameters"))
    return nothing
end

"""
    pulse_value(pulse::AbstractPulse, t, duration)

Evaluate a pulse at physical time `t` for an owning gate of nonnegative finite
`duration`. Parameterized pulses require a finite scalar result. Built-in pulses
return their offset outside their active window. `pulse(t, duration)` is equivalent.
Implement this method for custom `AbstractPulse` subtypes.
"""
function pulse_value end

"""Validate a pulse against its owning gate duration."""
validate_pulse(::AbstractPulse, duration) = nothing

(pulse::AbstractPulse)(t, duration) = pulse_value(pulse, t, duration)

function _require_finite(name, value)
    value isa Number && isfinite(value) || throw(ArgumentError("$name must be a finite number"))
    return value
end

function _require_symbol(name, value)
    value isa Symbol || throw(ArgumentError("$name must be a Symbol"))
    return value
end

function _require_finite_real(name, value)
    value isa Real && isfinite(value) || throw(ArgumentError("$name must be a finite real number"))
    return value
end

function pulse_value(pulse::AbstractParameterizedPulse, t, duration)
    _require_finite_real("duration", duration)
    duration >= 0 || throw(ArgumentError("Gate duration must be nonnegative"))
    parameters = merge(pulse.parameters, (; duration))
    return _require_finite("Pulse value", pulse_function(pulse)(parameters, t))
end

function validate_pulse(pulse::AbstractParameterizedPulse, duration)
    pulse_value(pulse, 0.0, duration)
    pulse_value(pulse, duration, duration)
    return nothing
end

include("internal_pulse_function.jl")
include("generic_pulse_function.jl")
include("pulses/constant_pulse.jl")
include("pulses/gaussian_pulse.jl")
include("pulses/sine_squared_pulse.jl")
include("pulses/sine_pulse.jl")
include("pulses/ramped_flattop.jl")
