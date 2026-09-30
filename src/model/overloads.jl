function Base.getproperty(model::DeviceModel, id::Symbol)
    hasfield(DeviceModel, id) && return getfield(model, id)
    hasfield(TrackingResult, id) && return getproperty(getfield(model, :eigensystem), id)
    throw(ArgumentError("Property $id not found in DeviceModel"))
end

Base.propertynames(::DeviceModel, private::Bool = false) =
    (fieldnames(DeviceModel)..., fieldnames(TrackingResult)...)
function _replace_property(model::DeviceModel, key::Symbol, value)
    if key === :gates
        fields = map(fieldnames(DeviceModel)) do name
            name === :gates ? value : getfield(model, name)
        end
        return DeviceModel(fields...)
    end
    key in (:components, :interactions, :coupling_parameters, :max_dimension, :truncation_dimensions) ||
        throw(ArgumentError("DeviceModel.$key is derived or unknown"))
    components = key === :components ? value : model.components
    length(components) == length(model.components) &&
        all(a.name == b.name for (a, b) in zip(components, model.components)) ||
        throw(ArgumentError("Component updates must preserve names and order"))
    truncations = if key === :truncation_dimensions
        value
    elseif key === :components
        Dict(new => model.truncation_dimensions[old] for (old, new) in zip(model.components, components)
            if haskey(model.truncation_dimensions, old))
    else
        model.truncation_dimensions
    end
    result = make_model(components,
        key === :interactions ? value : model.interactions,
        key === :coupling_parameters ? value : model.coupling_parameters;
        max_dimension = key === :max_dimension ? value : model.max_dimension,
        truncation_dimensions = truncations)
    merge!(result.gates, model.gates)
    return result
end
