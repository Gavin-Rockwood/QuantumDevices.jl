const _MODEL_INPUT_FIELDS = (:components, :interactions, :coupling_parameters,
    :max_dimension, :truncation_dimensions, :gates)

function _copy_model_fields(model, changes)
    return DeviceModel((get(changes,key,getfield(model,key)) for key in fieldnames(DeviceModel))...)
end

function _rebuild_model(model::DeviceModel, changes; threaded_components::Bool=true,
    reference_model=nothing,tracking_steps=nothing)
    all(key -> key in _MODEL_INPUT_FIELDS,keys(changes)) ||
        throw(ArgumentError("Only model physical inputs and gates are writable"))
    changes = Dict(key=>value for (key,value) in changes if key === :gates || value !== getfield(model,key))
    isempty(changes) && return model
    all(key -> key === :gates,keys(changes)) && return _copy_model_fields(model,changes)
    components = get(changes,:components,model.components)
    length(components)==length(model.components) &&
        all(a.name==b.name for (a,b) in zip(components,model.components)) ||
        throw(ArgumentError("Component updates must preserve names and order"))
    truncations = get(changes,:truncation_dimensions,
        Dict(new=>model.truncation_dimensions[old] for (old,new) in zip(model.components,components)
            if haskey(model.truncation_dimensions,old)))
    if reference_model === :initial
        compatible = try
            _reference_transport(model,components,_model_dimensions(components,truncations,get(changes,:max_dimension,model.max_dimension)))
            true
        catch error
            error isa Union{ArgumentError,DimensionMismatch} || rethrow()
            false
        end
        reference_model = compatible ? model : nothing
        tracking_steps = compatible ? tracking_steps : nothing
    end
    result = make_model(components,get(changes,:interactions,model.interactions),
        get(changes,:coupling_parameters,model.coupling_parameters);
        max_dimension=get(changes,:max_dimension,model.max_dimension),
        truncation_dimensions=truncations,threaded_components,reference_model,tracking_steps)
    return haskey(changes,:gates) ? _copy_model_fields(result,Dict(:gates=>copy(changes[:gates]))) : result
end
