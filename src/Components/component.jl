"""
    Component(spec, id, uuid = _new_uuid(), name = String(id), metadata = Dict(); ...)

Concrete component instance with resolved parameter catalog, available
operators, and symbolic Hamiltonian expression.
"""
struct Component{CT <: AbstractComponentSpec}
    spec::CT
    id::Symbol
    uuid::UUID
    name::String
    parameters::Dict{Symbol, DeviceParameter}
    operators::Set{Any}
    hamiltonian::Union{Nothing, Expr}
    metadata::Dict{Symbol, Any}
end

"""
    Component(spec, id; hamiltonian = default_hamiltonian(spec), extra_terms = nothing)
    Component(id, SpecType, args...; hamiltonian = missing, extra_terms = nothing)

Create a component from a spec, collecting parameters referenced by the stored
Hamiltonian.
"""
function Component(
    spec::AbstractComponentSpec,
    id::Symbol,
    uuid::UUID = _new_uuid(),
    name::AbstractString = String(id),
    metadata = Dict{Symbol, Any}();
    hamiltonian = default_hamiltonian(spec),
    extra_terms = nothing,
)
    component_hamiltonian = _component_hamiltonian(hamiltonian, extra_terms)
    
    return Component(
        spec,
        id,
        uuid,
        String(name),
        _component_parameters(spec, component_hamiltonian),
        Set{Any}(available_operators(spec)),
        component_hamiltonian,
        _metadata_dict(metadata),
    )
end

function Component(
    id,
    SpecType::Type{<:AbstractComponentSpec},
    args...;
    name = String(id),
    uuid = _new_uuid(),
    metadata = Dict{Symbol, Any}(),
    hamiltonian = missing,
    extra_terms = nothing,
)
    spec = SpecType(args...)
    component_hamiltonian =
        hamiltonian === missing ? default_hamiltonian(spec) : hamiltonian
    return Component(
        spec,
        id,
        uuid,
        name,
        _metadata_dict(metadata);
        hamiltonian = component_hamiltonian,
        extra_terms = extra_terms,
    )
end

"""
    default_hamiltonian(spec)

Return `spec.hamiltonian` when the specification defines one.
"""
default_hamiltonian(spec::AbstractComponentSpec) =
    hasproperty(spec, :hamiltonian) ? getproperty(spec, :hamiltonian) : nothing

"""
    _component_hamiltonian(hamiltonian, extra_terms)

Combine a component's base Hamiltonian expression with optional extra terms.
"""
function _component_hamiltonian(hamiltonian, extra_terms)
    if hamiltonian === nothing
        return extra_terms
    elseif extra_terms === nothing
        return hamiltonian
    else
        return hamiltonian + extra_terms
    end
end

"""
    _component_parameters(spec, expr)

Collect spec parameters plus any `DeviceParameter`s referenced by a Hamiltonian
expression.
"""
function _component_parameters(spec::AbstractComponentSpec, expr)
    parameters = _spec_parameters(spec)
    expr === nothing && return parameters

    for parameter in expression_parameters(expr)
        if parameter isa DeviceParameter
            parameters[_parameter_key(parameter)] = parameter
        elseif parameter isa Tuple
            key = _path_key(parameter)
            haskey(parameters, key) || error("Parameter $(join(parameter, ".")) is not defined on component spec $(typeof(spec)).")
        end
    end
    return parameters
end

"""
    _spec_parameters(spec)

Extract `DeviceParameter` fields declared directly on a component spec.
"""
function _spec_parameters(spec::AbstractComponentSpec)
    parameters = Dict{Symbol, DeviceParameter}()
    for field in fieldnames(typeof(spec))
        value = getfield(spec, field)
        value isa DeviceParameter || continue
        parameters[_parameter_key(value)] = value
    end
    return parameters
end

_parameter_key(parameter::DeviceParameter) = _path_key(parameter.path.parts)
_path_key(path::Tuple{Vararg{Symbol}}) = path[end]

"""
    expression_parameters(expr)

Return parameter references contained in a scalar or operator expression.
"""
function expression_parameters(expr::ScalarExpr, parameters = Set{Any}())
    if expr.head == :param
        push!(parameters, expr.args[1])
    else
        foreach(arg -> arg isa Expr && expression_parameters(arg, parameters), expr.args)
    end
    return parameters
end

function expression_parameters(expr::OperatorExpr, parameters = Set{Any}())
    foreach(arg -> arg isa Expr && expression_parameters(arg, parameters), expr.args)
    return parameters
end

"""
    numerical(component; dimension = nothing, kwargs...)
    numerical(component, operator; dimension = nothing, kwargs...)
    numerical(component, expr; dimension = nothing, kwargs...)

Evaluate a component Hamiltonian, catalog operator, or symbolic expression.
"""
function numerical(component::Component; dimension = nothing, kwargs...)
    component.hamiltonian === nothing && error("Component $(component.id) does not define a Hamiltonian.")
    return numerical(component, component.hamiltonian; dimension = _numerical_dimension(component, dimension), kwargs...)
end

function numerical(component::Component, operator::Symbol; dimension = nothing, kwargs...)
    operator in component.operators || error("Operator $operator is not available on component $(component.id).")
    return numerical_operator(component, operator; dimension = _numerical_dimension(component, dimension), kwargs...)
end

function numerical(component::Component, operator::Tuple{Vararg{Symbol}}; dimension = nothing, kwargs...)
    operator in component.operators || error("Operator $(join(operator, "/")) is not available on component $(component.id).")
    return numerical_operator(component, operator; dimension = _numerical_dimension(component, dimension), kwargs...)
end

"""
    numerical_operator(component, operator; dimension = nothing, kwargs...)

Forward operator evaluation from a component instance to its specification.
"""
numerical_operator(component::Component, operator; dimension = nothing, kwargs...) =
    numerical_operator(component.spec, operator; dimension = dimension, kwargs...)

function numerical(component::Component, expr::OperatorExpr; dimension = nothing, kwargs...)
    resolved_dimension = _numerical_dimension(component, dimension)
    if expr.head == :op
        operator = _lookup_operator(component, expr.args[1].path)
        return numerical(component, operator; dimension = resolved_dimension, kwargs...)
    elseif expr.head == :+
        return reduce(+, (numerical(component, arg; dimension = resolved_dimension, kwargs...) for arg in expr.args))
    elseif expr.head == :*
        return reduce(*, (_numerical_factor(component, arg; dimension = resolved_dimension, kwargs...) for arg in expr.args))
    elseif expr.head == :^
        return numerical(component, expr.args[1]; dimension = resolved_dimension, kwargs...) ^ expr.args[2]
    elseif expr.head == :adjoint
        return adjoint(numerical(component, expr.args[1]; dimension = resolved_dimension, kwargs...))
    else
        error("Unknown operator expression head $(expr.head).")
    end
end

function numerical(component::Component, expr::ScalarExpr; kwargs...)
    if expr.head == :const
        return expr.args[1]
    elseif expr.head == :param
        return _parameter_value(component, expr.args[1]; kwargs...)
    elseif expr.head == :+
        return reduce(+, (numerical(component, arg; kwargs...) for arg in expr.args))
    elseif expr.head == :*
        return reduce(*, (numerical(component, arg; kwargs...) for arg in expr.args))
    elseif expr.head == :/
        return numerical(component, expr.args[1]; kwargs...) / numerical(component, expr.args[2]; kwargs...)
    elseif expr.head == :^
        return numerical(component, expr.args[1]; kwargs...) ^ expr.args[2]
    elseif expr.head == :neg
        return -numerical(component, expr.args[1]; kwargs...)
    elseif expr.head == :cos
        return cos(numerical(component, expr.args[1]; kwargs...))
    elseif expr.head == :sin
        return sin(numerical(component, expr.args[1]; kwargs...))
    else
        error("Unknown scalar expression head $(expr.head).")
    end
end

_numerical_factor(component::Component, arg::OperatorExpr; dimension, kwargs...) =
    numerical(component, arg; dimension = dimension, kwargs...)

_numerical_factor(component::Component, arg::ScalarExpr; dimension = nothing, kwargs...) =
    numerical(component, arg; kwargs...)

_numerical_factor(::Component, arg::Number; kwargs...) = arg

"""
    _numerical_dimension(component, dimension)

Use an explicit numerical dimension or the component's finite default dimension.
"""
function _numerical_dimension(component::Component, dimension)
    dimension !== nothing && return dimension
    haskey(component.parameters, :dimension) ||
        error("Component $(component.id) does not define a dimension parameter.")
    default = component.parameters[:dimension].default
    isfinite(default) ||
        error("Component $(component.id) has infinite maximum dimension; provide dimension explicitly.")
    return Int(default)
end

"""
    _lookup_operator(component, path)

Resolve an operator expression path to a key in the component's operator set.
"""
function _lookup_operator(component::Component, path::Tuple{Vararg{Symbol}})
    key = _relative_catalog_key(component, path, :operators)
    key in component.operators || error("Operator $key is not available on component $(component.id).")
    return key
end

_parameter_key(parameter::DeviceParameter, component::Component) =
    _relative_catalog_key(component, parameter.path.parts, :params)

_parameter_key(path::Tuple{Vararg{Symbol}}, component::Component) =
    _relative_catalog_key(component, path, :params)

"""
    _parameter_value(component, parameter; kwargs...)

Resolve a component parameter from keyword overrides, `params`, or defaults.
"""
function _parameter_value(component::Component, parameter; kwargs...)
    key = _parameter_key(parameter, component)
    values = Dict{Symbol, Any}(pairs(kwargs))
    if haskey(values, key)
        return values[key]
    end

    if haskey(values, :params)
        params = values[:params]
        if haskey(params, key)
            return params[key]
        elseif parameter isa DeviceParameter && haskey(params, parameter.path)
            return params[parameter.path]
        end
    end

    parameter isa DeviceParameter && return parameter.default
    haskey(component.parameters, key) ||
        error("Parameter $key is not available on component $(component.id).")
    return component.parameters[key].default
end

"""
    _relative_catalog_key(component, path, catalog)

Convert absolute or relative expression paths into component-local catalog keys.
"""
function _relative_catalog_key(component::Component, path::Tuple{Vararg{Symbol}}, catalog::Symbol)
    if length(path) >= 4 && path[1] == :components
        path[2] == component.id || error("Path $(join(path, "/")) does not belong to component $(component.id).")
        path[3] == catalog || error("Path $(join(path, "/")) does not refer to $catalog.")
        relative = path[4:end]
        return length(relative) == 1 ? relative[1] : relative
    elseif length(path) >= 2 && path[1] == catalog
        relative = path[2:end]
        return length(relative) == 1 ? relative[1] : relative
    else
        return length(path) == 1 ? path[1] : path
    end
end
