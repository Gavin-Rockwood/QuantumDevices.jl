"""
    numerical(model::DeviceModel, gate::DeviceGate)

Combine the projected idle Hamiltonian and gate drive with gate parameter overrides.
Overrides of idle parameters must equal their idle values at both endpoints.
"""
function numerical(model::DeviceModel, gate::DeviceGate)
    combined = Dict{Symbol,Any}(pairs(model.parameters))
    for (key, value) in pairs(gate.parameters)
        if haskey(combined, key) && !_check_gate_endpoints(value, combined[key], gate.duration)
            throw(ArgumentError("Endpoints of parameter $key do not agree with the idle parameter"))
        end
        combined[key] = value isa AbstractPulse ?
            (t -> _pulse_scalar_value(value, t, gate.duration)) : value
    end
    required = union(parameter_keys(model.hamiltonian), parameter_keys(gate.hamiltonian))
    missing_ = setdiff(required, Set(keys(combined)))
    isempty(missing_) || throw(ArgumentError("Missing gate parameters: $missing_"))
    dims = [get(model.truncation_dimensions, c, c.dimension) for c in model.components]
    drive_terms = _projected_terms(gate.hamiltonian, model.components, dims)
    return _evaluate_terms(vcat(model.compiled_terms, drive_terms), combined)
end
