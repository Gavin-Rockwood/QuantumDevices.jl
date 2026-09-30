"""
    save(path::AbstractString, model::DeviceModel)

Save a model and its gates to a new directory. Component and interaction values
live in `parameters.json`; gate records live in `gates/`. Existing paths are never
replaced. Generic pulse artifacts require their callable definitions when loaded.
"""
function save(path::AbstractString, model::DeviceModel)
    destination = abspath(path)
    ispath(destination) && throw(ArgumentError("Destination already exists: $destination"))
    mkpath(dirname(destination))
    staging = mktempdir(dirname(destination); prefix = ".quantumdevices-")
    try
        _save_model(staging, model)
        mv(staging, destination)
    finally
        isdir(staging) && rm(staging; recursive = true)
    end
    return destination
end

function _save_model(directory, model::DeviceModel)
    component_parameters = Dict(c.name => _encode_parameters(c.parameters) for c in model.components)
    _write_json(joinpath(directory, "parameters.json"), Dict("schema_version" => 1,
        "type" => "parameters", "components" => component_parameters,
        "interactions" => _encode_parameters(model.coupling_parameters)))
    components = [_component_record(c, joinpath(directory, "components"), i)
        for (i, c) in enumerate(model.components)]
    gates = map(collect(enumerate(pairs(model.gates)))) do (index, (key, gate))
        gate isa DeviceGate || throw(ArgumentError("Model gates must contain DeviceGate values"))
        folder = lpad(string(index), 4, '0')
        _save_gate(joinpath(directory, "gates", folder), gate)
        Dict("key" => _encode_value(key), "directory" => folder)
    end
    _write_json(joinpath(directory, "model.json"), Dict("schema_version" => 1, "type" => "model",
        "components" => components, "interactions" => _expression_record(model.interactions),
        "max_dimension" => model.max_dimension, "gates" => gates,
        "truncation_dimensions" => Dict(c.name => d for (c, d) in model.truncation_dimensions)))
end

"""
    load(path::AbstractString)

Load the model bundle directory at `path` and return a [`DeviceModel`](@ref).
Reconstruct component bases, promoted operators, the idle Hamiltonian, and dressed
states; then restore named gates, preserving string versus symbol gate keys.
Requires `model.json`, `parameters.json`, and all referenced artifacts.

Unsupported schema tags, stored values, or symbolic operations throw. Generic
pulse callable definitions must be available in the process before loading.
Use trusted bundles, especially those containing JLD2 artifacts. This method
loads the directory format created by [`save`](@ref), not legacy source archives.
"""
function load(path::AbstractString)
    data = _read_json(joinpath(path, "model.json"))
    _check_schema(data, "model")
    parameters = _read_json(joinpath(path, "parameters.json"))
    _check_schema(parameters, "parameters")
    components = [_restore_component(c, parameters["components"][c["name"]], joinpath(path, "components"))
        for c in data["components"]]
    truncations = Dict(c => data["truncation_dimensions"][c.name] for c in components
        if haskey(data["truncation_dimensions"], c.name))
    model = make_model(components, _restore_expression(data["interactions"]),
        _decode_parameters(parameters["interactions"]); max_dimension = data["max_dimension"],
        truncation_dimensions = truncations)
    for record in data["gates"]
        key = _decode_value(record["key"])
        haskey(model.gates, key) && throw(ArgumentError("Duplicate stored gate key: $key"))
        model.gates[key] = _load_gate(_bundle_child(joinpath(path, "gates"), record["directory"]))
    end
    return model
end
