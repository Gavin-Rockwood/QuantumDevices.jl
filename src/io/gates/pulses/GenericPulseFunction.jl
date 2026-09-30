function _pulse_record(pulse::GenericPulseFunction, directory, index)
    mkpath(directory)
    filename = string(index, ".jld2")
    JLD2.save_object(joinpath(directory, filename), pulse.callable)
    return Dict("type" => "generic", "signature" => "parameters_time",
        "parameters" => _encode_parameters(pulse.parameters), "artifact" => filename)
end

function _restore_pulse(::Val{:generic}, data, directory)
    get(data, "signature", nothing) == "parameters_time" || throw(ArgumentError(
        "Unsupported generic pulse format; expected an f(p, t) callable with named parameters"))
    parameters = _decode_parameters(data["parameters"])
    callable = JLD2.load_object(_bundle_child(directory, data["artifact"]))
    applicable(callable, merge(parameters, (; duration = 1.0)), 0.0) || throw(ArgumentError(
        "Saved pulse callable could not be restored. Load its defining module before loading the bundle; arbitrary interactive functions are not portable."))
    return GenericPulseFunction(callable, parameters)
end
