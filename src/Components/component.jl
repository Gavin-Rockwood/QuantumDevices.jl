"""A named instance of a component specification."""
struct Component{S<:AbstractComponentSpec}
    spec::S
    name::Symbol
end

"""
    copy(component; name=component.name, params=nothing)

Create an independent component whose parameter defaults may be overridden.
`params` accepts mapping-like objects keyed by local parameter paths.
"""
function Base.copy(component::Component; name = component.name, params = nothing)
    name isa Symbol || throw(ArgumentError("Copied component name must be a Symbol."))
    overrides = _component_parameter_overrides(component, params)
    replacements = _component_parameter_copies(component, overrides)
    Component(_copy_component_spec(component.spec, replacements), name)
end

"""
    update!(component; params=nothing)

Update a component's parameter defaults in place. All supplied values are
validated before any default is changed.
"""
function update!(component::Component; params = nothing)
    overrides = _component_parameter_overrides(component, params)
    updates = Tuple{DeviceParameter,Any}[
        (parameter, overrides[parameter.path])
        for parameter in _component_parameter_objects(component)
        if haskey(overrides, parameter.path)
    ]
    for (parameter, value) in updates
        _parameter_value(parameter, value)
    end
    for (parameter, value) in updates
        parameter.default = value
    end
    component
end

function _component_parameter_overrides(component, params)
    overrides = Dict{ParamPath,Any}()
    params === nothing && return overrides

    supplied_keys = try
        keys(params)
    catch exception
        exception isa MethodError || rethrow()
        error("params must support keys(params) and params[key].")
    end

    available = parameters(component)
    for key in supplied_keys
        path = try
            _param_path(key)
        catch exception
            exception isa MethodError || rethrow()
            error("Unsupported component parameter key $key.")
        end
        haskey(available, path) || error(
            "Component $(component.name) does not define parameter $(join(path.parts, "/")).",
        )
        haskey(overrides, path) && error(
            "Component parameter path $(join(path.parts, "/")) was supplied more than once.",
        )
        value = try
            params[key]
        catch exception
            exception isa MethodError || rethrow()
            error("params must support keys(params) and params[key].")
        end
        _parameter_value(available[path], value)
        overrides[path] = value
    end
    overrides
end

function _component_parameter_objects(component)
    originals = DeviceParameter[]
    expression = hamiltonian(component)
    if expression isa Expr
        for parameter in get_params(expression)
            any(candidate -> candidate === parameter, originals) ||
                push!(originals, parameter)
        end
    end
    for field in fieldnames(typeof(component.spec))
        value = getfield(component.spec, field)
        value isa DeviceParameter || continue
        any(candidate -> candidate === value, originals) || push!(originals, value)
    end
    originals
end

function _component_parameter_copies(component, overrides)
    replacements = IdDict{DeviceParameter,DeviceParameter}()
    for parameter in _component_parameter_objects(component)
        default = get(overrides, parameter.path, parameter.default)
        replacements[parameter] = DeviceParameter(
            parameter.path,
            parameter.domain,
            parameter.fixed,
            parameter.required,
            default,
            parameter.description,
            copy(parameter.metadata),
        )
    end
    replacements
end

_replace_component_parameters(value, replacements) = value
_replace_component_parameters(parameter::DeviceParameter, replacements) =
    get(replacements, parameter, parameter)
function _replace_component_parameters(expression::Expr, replacements)
    args = map(
        argument -> _replace_component_parameters(argument, replacements),
        expression.args,
    )
    typeof(expression)(expression.head, args)
end

function _copy_component_spec(spec::AbstractComponentSpec, replacements)
    values = ntuple(fieldcount(typeof(spec))) do index
        _replace_component_parameters(getfield(spec, index), replacements)
    end
    try
        typeof(spec)(values...)
    catch exception
        exception isa MethodError || rethrow()
        error(
            "Cannot copy component spec $(nameof(typeof(spec))): " *
            "its concrete type must support construction from its fields.",
        )
    end
end

"""Return the symbolic Hamiltonian owned by a component specification."""
hamiltonian(spec::AbstractComponentSpec) = spec.hamiltonian
hamiltonian(component::Component) = hamiltonian(component.spec)

"""Return the operator names exposed by a component specification."""
operators(spec::AbstractComponentSpec) = Set{Any}(keys(spec.operators))
operators(component::Component) = operators(component.spec)

"""Collect the parameters carried by an expression."""
function parameters(expression::Expr)
    registry = Dict{ParamPath,DeviceParameter}()
    for parameter in get_params(expression)
        registry[parameter.path] = parameter
    end
    registry
end

"""Collect the parameters carried by a specification's Hamiltonian."""
function parameters(spec::AbstractComponentSpec)
    expression = hamiltonian(spec)
    expression isa Expr ? parameters(expression) : Dict{ParamPath,DeviceParameter}()
end

parameters(component::Component) = parameters(component.spec)

"""Evaluate a component's Hamiltonian."""
function numerical(component::Component; dimension = nothing, kwargs...)
    numerical(component, hamiltonian(component);
        dimension = _numerical_dimension(component, dimension), kwargs...)
end

function numerical(component::Component, operator::Union{Symbol,Tuple{Vararg{Symbol}}};
        dimension = nothing, kwargs...)
    dimension = _numerical_dimension(component, dimension)
    key = _operator_key(operator)
    key in operators(component) ||
        error("Operator $key is not available on component $(component.name).")
    numerical_operator(component.spec, key; dimension, kwargs...)
end

function numerical(component::Component, expression::OperatorExpr;
        dimension = nothing, kwargs...)
    dimension = _numerical_dimension(component, dimension)
    _numerical_expression(
        expression,
        path -> numerical(component, _operator_key(path); dimension, kwargs...),
        parameter -> _component_value(parameter; kwargs...),
    )
end

numerical(component::Component, expression::ScalarExpr; kwargs...) =
    _numerical_scalar(expression, parameter -> _component_value(parameter; kwargs...))

function _numerical_dimension(component::Component, dimension)
    dimension !== nothing && return _dimension_tuple(dimension)
    default = component.spec.dimension.size
    all(isfinite, default) ||
        error("Component $(component.name) has infinite maximum dimension; provide dimension explicitly.")
    default
end

_dimension_tuple(dimension::Tuple) = dimension
_dimension_tuple(dimension) =
    throw(ArgumentError("numerical dimension must be a tuple; got $(typeof(dimension))."))

_operator_key(operator::Symbol) = operator
_operator_key(operator::Tuple{Vararg{Symbol}}) = begin
    relative = _canonical_reference(operator, :operators)
    length(relative) == 1 ? only(relative) : relative
end

function _component_parameter(component::Component, name::Symbol)
    matches = [parameter for parameter in values(parameters(component))
        if parameter.path.parts == (name,)]
    isempty(matches) && error("Component $(component.name) does not define parameter $name.")
    length(matches) == 1 || error("Component $(component.name) has multiple parameters named $name.")
    only(matches)
end

function _component_value(parameter; kwargs...)
    parameter isa DeviceParameter ||
        error("Scalar expressions must carry DeviceParameter objects.")
    key = length(parameter.path.parts) == 1 ? only(parameter.path.parts) : parameter.path.parts
    values = Dict{Symbol,Any}(pairs(kwargs))
    key isa Symbol && haskey(values, key) && return _parameter_value(parameter, values[key])
    if haskey(values, :params)
        params = values[:params]
        haskey(params, key) && return _parameter_value(parameter, params[key])
        haskey(params, parameter.path) && return _parameter_value(parameter, params[parameter.path])
    end
    _parameter_value(parameter)
end
