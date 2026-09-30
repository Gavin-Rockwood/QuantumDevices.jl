function _component_record(component::Component, artifact_dir, index)
    data = Dict{String,Any}("name" => component.name, "type" => component.type,
        "dimension" => component.dimension)
    if !haskey(COMPONENT_CONSTRUCTORS, component.type)
        mkpath(artifact_dir)
        filename = string(index, ".jld2")
        JLD2.save_object(joinpath(artifact_dir, filename), component)
        data["artifact"] = filename
    end
    return data
end

function _restore_component(data, parameters, artifact_dir)
    name, kind, dimension = String(data["name"]), data["type"], data["dimension"]
    p = _decode_parameters(parameters)
    if haskey(data, "artifact")
        component = JLD2.load_object(_bundle_child(artifact_dir, data["artifact"]))
        component isa Component || throw(ArgumentError("Invalid component artifact"))
        return Component(name, p, component.operators, component.hamiltonian, String(kind), dimension)
    end
    constructor = get(COMPONENT_CONSTRUCTORS, String(kind), nothing)
    constructor === nothing && throw(ArgumentError("Unknown component type: $kind"))
    return constructor(; name, dimension, p...)
end
