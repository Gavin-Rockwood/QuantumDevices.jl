"""A named instance of a component specification."""
struct Component{S<:AbstractComponentSpec}
    spec::S
    name::Symbol
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
