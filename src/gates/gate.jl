"""
    DeviceGate(parameters::NamedTuple, hamiltonian, duration=nothing)

A gate with static numbers or timed [`Pulse`](@ref) controls. Its symbolic
Hamiltonian adds to the model's idle Hamiltonian. Controls with existing model
parameter names replace those values and must match idle values at both endpoints.

By default duration is the latest control end time and is recomputed after pulse
updates. An explicit finite nonnegative duration must contain every pulse and can
add trailing idle time. Scalar-only gates require an explicit duration.
Coefficients use frequency units; evolution requires `2π * H`.
"""
struct DeviceGate{P<:NamedTuple,H,D}
    parameters::P
    hamiltonian::H
    duration_override::D
    function DeviceGate(parameters::P, hamiltonian::H, duration::D) where {P<:NamedTuple,H,D}
        if duration !== nothing
            _require_finite_real("Gate duration", duration) >= 0 || throw(ArgumentError("Gate duration must be nonnegative"))
        end
        pulses = AbstractPulse[]
        for (name, value) in pairs(parameters)
            value isa Union{Number,AbstractPulse} || throw(ArgumentError("Gate parameter $name must be a number or timed pulse"))
            value isa Number && _require_finite("Gate parameter $name", value)
            if value isa AbstractPulse
                end_time = _validate_timed_control(value)
                push!(pulses, value)
                duration !== nothing && end_time > duration &&
                    throw(ArgumentError("Pulse $name ends beyond the explicit gate duration"))
            end
        end
        duration === nothing && isempty(pulses) && throw(ArgumentError("Scalar-only gates require an explicit duration"))
        new{P,H,D}(parameters, hamiltonian, duration)
    end
end
DeviceGate(parameters::NamedTuple, hamiltonian) = DeviceGate(parameters, hamiltonian, nothing)
