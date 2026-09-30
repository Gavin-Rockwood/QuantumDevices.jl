"""
    DeviceGate(parameters::NamedTuple, hamiltonian, duration::Real)

A gate with scalar or [`AbstractPulse`](@ref) parameters, an additional symbolic
Hamiltonian (or scalar zero), and finite nonnegative duration. Pulses are validated
at the endpoints. Evaluate it in a model with [`numerical`](@ref); its Hamiltonian
adds to the idle Hamiltonian, and parameters override model parameters.

Overrides of idle parameters must match their idle values exactly at both gate
endpoints. Raw functions are rejected: wrap them in [`GenericPulseFunction`](@ref).
Use [`setpath`](@ref) for immutable updates.
"""
struct DeviceGate{P<:NamedTuple,H,D<:Real}
    parameters::P
    hamiltonian::H
    duration::D

    function DeviceGate(parameters::P, hamiltonian::H, duration::D) where {P<:NamedTuple,H,D<:Real}
        isfinite(duration) && duration >= 0 ||
            throw(ArgumentError("Gate duration must be finite and nonnegative"))
        for (name, value) in pairs(parameters)
            value isa Union{Number,AbstractPulse} ||
                throw(ArgumentError("Gate parameter $name must be a Number or AbstractPulse"))
            value isa Number && !isfinite(value) &&
                throw(ArgumentError("Gate parameter $name must be finite"))
            value isa AbstractPulse && validate_pulse(value, duration)
        end
        new{P,H,D}(parameters, hamiltonian, duration)
    end
end

function _pulse_scalar_value(pulse::AbstractPulse, t, duration)
    value = pulse_value(pulse, t, duration)
    value isa Number || throw(ArgumentError("Pulse functions must return scalar numbers"))
    return value
end

function _check_gate_endpoints(gate_param, model_param, duration)
    return gate_param == model_param
end
function _check_gate_endpoints(gate_param::AbstractPulse, model_param, duration)
    return (_pulse_scalar_value(gate_param, 0.0, duration) == model_param) &&
           (_pulse_scalar_value(gate_param, duration, duration) == model_param)
end
