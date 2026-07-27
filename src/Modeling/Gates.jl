"""A symbolic, finite-duration Hamiltonian overlay for one exact ModelSpec."""
struct GateSpec{M<:ModelSpec}
    name::Symbol
    modelspec::M
    duration::Float64
    hamiltonian::Union{Nothing,OperatorExpr}
    parameters::Dict{ParamPath,DeviceParameter}
end

function GateSpec(
        name::Symbol,
        modelspec::ModelSpec;
        duration::Real,
        hamiltonian = nothing,
        parameters = ())
    duration > 0 || error("GateSpec duration must be positive.")
    hamiltonian === nothing || hamiltonian isa OperatorExpr ||
        error("GateSpec hamiltonian must be an OperatorExpr or nothing.")

    registry = Dict{ParamPath,DeviceParameter}()
    for parameter in parameters
        parameter isa DeviceParameter ||
            error("GateSpec parameters must be DeviceParameter objects.")
        registry[parameter.path] = parameter
    end
    for (path, parameter) in registry
        haskey(modelspec.parameters, path) || continue
        modelspec.parameters[path].fixed &&
            error("GateSpec $name cannot override fixed model parameter $(join(path.parts, "/")).")
    end
    if hamiltonian !== nothing
        for (path, parameter) in QuantumDevices.parameters(hamiltonian)
            haskey(registry, path) && continue
            try
                _model_parameter_path(modelspec, parameter)
            catch
                error("Parameter $(join(path.parts, "/")) in gate $name is not registered by the gate or its ModelSpec.")
            end
        end
    end
    GateSpec(name, modelspec, Float64(duration), hamiltonian, registry)
end

"""Build the complete static or time-dependent Hamiltonian for a gate."""
function numerical(model::QuantumDeviceModel, gate::GateSpec)
    model.spec === gate.modelspec ||
        error("GateSpec $(gate.name) belongs to a different ModelSpec.")

    active = _gate_model_hamiltonian(model, gate, false)
    baseline = _gate_model_hamiltonian(model, gate, true)
    result = model.hamiltonian + active - baseline
    if gate.hamiltonian !== nothing
        result += _numerical_expression(
            gate.hamiltonian,
            path -> _gate_operator(model, path),
            parameter -> _gate_expression_value(model.spec, gate, parameter),
        )
    end

    _check_gate_endpoint(model, gate, result, 0.0)
    _check_gate_endpoint(model, gate, result, gate.duration)
    result
end

numerical(gate::GateSpec) = numerical(model(gate.modelspec), gate)

function _gate_model_hamiltonian(model, gate, baseline)
    total = nothing
    for (name, component) in pairs(model.spec.subsystems)
        value = _numerical_expression(
            hamiltonian(component),
            path -> _gate_operator(model, path, name),
            parameter -> begin
                path = ParamPath(name, parameter.path.parts...)
                _gate_model_value(model.spec, gate, path, baseline)
            end,
        )
        total = total === nothing ? value : total + value
    end
    for interaction in model.spec.interactions
        value = _numerical_expression(
            interaction,
            path -> _gate_operator(model, path),
            parameter -> _gate_model_value(
                model.spec,
                gate,
                _model_parameter_path(model.spec, parameter),
                baseline,
            ),
        )
        total = total === nothing ? value : total + value
    end
    total
end

function _gate_model_value(spec, gate, path, baseline)
    value = if haskey(gate.parameters, path)
        _parameter_value(spec.parameters[path], _parameter_value(gate.parameters[path]))
    else
        _model_value(spec, path)
    end
    baseline ? _at(_model_value(spec, path), 0) : value
end

function _gate_expression_value(spec, gate, parameter)
    haskey(gate.parameters, parameter.path) &&
        return _parameter_value(gate.parameters[parameter.path])
    path = _model_parameter_path(spec, parameter)
    haskey(gate.parameters, path) &&
        return _parameter_value(spec.parameters[path], _parameter_value(gate.parameters[path]))
    _model_value(spec, path)
end

function _gate_operator(model, path, component = nothing)
    rooted = length(path) >= 1 && path[1] == :components
    relative = _canonical_reference(path, :operators)
    key = rooted || component === nothing ? relative : (component, relative...)
    haskey(model.operators, key) ||
        error("Operator path $(join(key, "/")) is not available in model $(model.spec.name).")
    model.operators[key]
end

function _check_gate_endpoint(model, gate, hamiltonian, time)
    value = hamiltonian isa QuantumObjectEvolution ? hamiltonian(nothing, time) : hamiltonian
    residual = norm(value - model.hamiltonian)
    tolerance = 1e-10 + 1e-10 * norm(model.hamiltonian)
    residual <= tolerance || error(
        "GateSpec $(gate.name) does not return to the model Hamiltonian at t=$time " *
        "(residual $residual, tolerance $tolerance).",
    )
end
