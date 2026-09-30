function _pulse_record(pulse::InternalPulseFunction, directory, index)
    return Dict("type" => "internal", "name" => String(pulse.name),
        "parameters" => _encode_parameters(pulse.parameters))
end

_restore_pulse(::Val{:internal}, data, directory) =
    InternalPulseFunction(Symbol(data["name"]), _decode_parameters(data["parameters"]))
