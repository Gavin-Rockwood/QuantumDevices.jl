"""
    numerical(model::DeviceModel, expression::Sym; scalar=1, dense=false)

Evaluate a symbolic expression using the model's component operators, retained
dimensions, and merged parameters. Operator names use the model's promoted names
(for example, `op(:q_charge)`). Complete operator products are multiplied before
projection, following the same rules as [`make_model`](@ref).

For example, `numerical(model, model.hamiltonian)` reproduces `model.H`.
Unresolved parameters produce a `QobjEvo`, as with the other symbolic overloads.
`dense=true` materializes all operator matrices as dense matrices at construction.
`scalar` is folded into operator matrices and constant coefficients at construction.
"""
function numerical(model::DeviceModel, expression::Sym; scalar::Number=1, dense::Bool=false)
    dims = [get(model.truncation_dimensions, c, c.dimension) for c in model.components]
    return _evaluate_terms(_projected_terms(expression, model.components, dims), model.parameters; scalar, dense)
end

function Base.getproperty(model::DeviceModel, id::Symbol)
    hasfield(DeviceModel, id) && return getfield(model, id)
    id === :spectrum && return getfield(model, :eigensystem).others
    hasfield(TrackingResult, id) && return getproperty(getfield(model, :eigensystem), id)
    throw(ArgumentError("Property $id not found in DeviceModel"))
end

Base.propertynames(::DeviceModel, private::Bool = false) =
    (fieldnames(DeviceModel)..., :spectrum, fieldnames(TrackingResult)...)
function _replace_property(model::DeviceModel, key::Symbol, value)
    return _rebuild_model(model, Dict(key => value))
end

_replace_properties(model::DeviceModel, changes) = _rebuild_model(model, changes)

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
