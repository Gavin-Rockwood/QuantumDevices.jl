"""
    numerical(model::DeviceModel, expression::Sym)

Evaluate a symbolic expression using the model's component operators, retained
dimensions, and merged parameters. Operator names use the model's promoted names
(for example, `op(:q_charge)`). Complete operator products are multiplied before
projection, following the same rules as [`make_model`](@ref).

For example, `numerical(model, model.hamiltonian)` reproduces `model.H`.
Unresolved parameters produce a `QobjEvo`, as with the other symbolic overloads.
"""
function numerical(model::DeviceModel, expression::Sym)
    dims = [get(model.truncation_dimensions, c, c.dimension) for c in model.components]
    return _evaluate_terms(_projected_terms(expression, model.components, dims), model.parameters)
end

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

function Base.show(io::IO, model::DeviceModel)
    print(io, "DeviceModel(", length(model.components), " components, dimension=",
        size(model.H, 1), ")")
end

function Base.show(io::IO, ::MIME"text/plain", model::DeviceModel)
    show(io, model)
    limited = IOContext(io, :compact => true, :limit => true)
    print(io, "\n  Components: ")
    show(limited, [c.name for c in model.components])
    print(io, "\n  Parameters: ")
    show(limited, model.parameters)
    print(io, "\n  Operators: ", length(model.operators))
    print(io, "\n  Compiled terms: ", length(model.compiled_terms))
    print(io, "\n  Gates: ")
    show(limited, collect(keys(model.gates)))
end
