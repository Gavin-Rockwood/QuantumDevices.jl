function _save_gate(directory, gate::DeviceGate)
    mkpath(directory)
    parameters = Dict{String,Any}()
    for (index, (name, value)) in enumerate(pairs(gate.parameters))
        parameters[String(name)] = value isa AbstractPulse ?
            Dict("kind" => "pulse", "pulse" => _control_record(value, joinpath(directory, "pulses"), index)) :
            Dict("kind" => "value", "value" => _encode_value(value))
    end
    _write_json(joinpath(directory, "gate.json"), Dict("schema_version" => 2, "type" => "gate",
        "duration" => gate.duration_override, "hamiltonian" => _expression_record(gate.hamiltonian),
        "parameters" => parameters))
end

function _load_gate(directory)
    data = _read_json(joinpath(directory, "gate.json"))
    data["type"] == "gate" || throw(ArgumentError("Expected gate record"))
    data["schema_version"] == 2 || throw(ArgumentError(
        "Unsupported legacy gate/pulse format; recreate gates with the timed Pulse API"))
    parameters = map(collect(pairs(data["parameters"]))) do (name, record)
        kind = record["kind"]
        value = if kind == "pulse"
            _restore_control(record["pulse"], joinpath(directory, "pulses"))
        elseif kind == "value"
            _decode_value(record["value"])
        else
            throw(ArgumentError("Unknown gate parameter kind: $kind"))
        end
        Symbol(name) => value
    end
    return DeviceGate((; parameters...), _restore_expression(data["hamiltonian"]), data["duration"])
end
