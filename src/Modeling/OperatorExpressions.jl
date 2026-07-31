"""
    OperatorRef(path)

Symbolic path to an operator in either a component catalog or a local component
catalog.
"""
struct OperatorRef
    path::Tuple{Vararg{Symbol}}
end

abstract type ExprRole end
struct ScalarRole{Data} <: ExprRole end
struct OperatorRole <: ExprRole end

"""
    Expr{R}

Minimal symbolic expression node used for scalar parameters and operator
algebra.
"""
struct Expr{R<:ExprRole}
    head::Symbol
    args::Tuple
end

const ScalarExpr = Expr{<:ScalarRole}
const OperatorExpr = Expr{OperatorRole}

ScalarExpr(head::Symbol, args::Tuple) = Expr{ScalarRole{Any}}(head, args)

constant(value::Number) = Expr{ScalarRole{typeof(value)}}(:const, (value,))

"""
    param(parameter)

Create a scalar parameter leaf that carries a `DeviceParameter` directly.
"""
param(parameter::DeviceParameter) =
    Expr{ScalarRole{DeviceParameter}}(:param, (parameter,))

"""
    op(name)
    op(component, path...)

Create an operator reference expression, either relative to a component or
rooted at another component in a model.
"""
op(name::Symbol) = OperatorExpr(:op, (OperatorRef((:operators, name)),))

function op(component::Symbol, path::Symbol...)
    isempty(path) && return op(component)
    return OperatorExpr(
        :op,
        (OperatorRef((:components, component, :operators, path...)),),
    )
end

_scalar(expr::ScalarExpr) = expr
_scalar(value::Number) = constant(value)
_scalar(value::DeviceParameter) = param(value)

import Base: +, -, *, /, ^, adjoint, conj, cos, sin

function +(a::ScalarExpr, b::Union{ScalarExpr,Number})
    return ScalarExpr(:+, (_flatten(:+, a)..., _flatten(:+, _scalar(b))...))
end

function +(a::Number, b::ScalarExpr)
    return _scalar(a) + b
end

function -(a::ScalarExpr, b::Union{ScalarExpr,Number})
    return ScalarExpr(:+, (_flatten(:+, a)..., -_scalar(b)))
end

function -(a::Number, b::ScalarExpr)
    return _scalar(a) - b
end

-(a::ScalarExpr) = ScalarExpr(:neg, (a,))

function *(a::ScalarExpr, b::Union{ScalarExpr,Number})
    return ScalarExpr(:*, (_flatten(:*, a)..., _flatten(:*, _scalar(b))...))
end

function *(a::Number, b::ScalarExpr)
    return _scalar(a) * b
end

/(a::ScalarExpr, b::Union{ScalarExpr,Number}) = ScalarExpr(:/, (a, _scalar(b)))
/(a::Number, b::ScalarExpr) = ScalarExpr(:/, (_scalar(a), b))
^(a::ScalarExpr, exponent::Number) = ScalarExpr(:^, (a, exponent))
cos(a::ScalarExpr) = ScalarExpr(:cos, (a,))
sin(a::ScalarExpr) = ScalarExpr(:sin, (a,))
conj(a::ScalarExpr) = ScalarExpr(:conj, (a,))

function +(a::OperatorExpr, b::OperatorExpr)
    return OperatorExpr(:+, (_flatten(:+, a)..., _flatten(:+, b)...))
end

function -(a::OperatorExpr, b::OperatorExpr)
    return OperatorExpr(:+, (_flatten(:+, a)..., -1 * b))
end

function *(a::OperatorExpr, b::OperatorExpr)
    return OperatorExpr(:*, (_flatten(:*, a)..., _flatten(:*, b)...))
end

function *(a::Union{ScalarExpr,Number}, b::OperatorExpr)
    return OperatorExpr(:*, (_scalar(a), _flatten(:*, b)...))
end

function *(a::OperatorExpr, b::Union{ScalarExpr,Number})
    return OperatorExpr(:*, (_flatten(:*, a)..., _scalar(b)))
end

function /(a::OperatorExpr, b::Union{ScalarExpr,Number})
    return OperatorExpr(:*, (_flatten(:*, a)..., constant(1) / _scalar(b)))
end

^(a::OperatorExpr, exponent::Number) = OperatorExpr(:^, (a, exponent))
adjoint(a::OperatorExpr) = OperatorExpr(:adjoint, (a,))

_flatten(head::Symbol, expr::Expr) = expr.head == head ? expr.args : (expr,)

"""
    get_params(expr)

Return the unique parameter references in a scalar or operator expression, in
first-occurrence order.
"""
function get_params(expr::Expr)
    parameters = DeviceParameter[]
    _get_params!(parameters, expr)
    return parameters
end

function _get_params!(parameters, expr::Expr)
    foreach(arg -> arg isa Expr && _get_params!(parameters, arg), expr.args)
    return parameters
end

function _get_params!(parameters, expr::Expr{ScalarRole{DeviceParameter}})
    parameter = only(expr.args)
    any(existing -> existing === parameter, parameters) || push!(parameters, parameter)
    return parameters
end

_get_params!(parameters, ::Expr{ScalarRole{Data}}) where {Data<:Number} = parameters

struct _ParameterizedCoefficient{F}
    evaluate::F
end

(coefficient::_ParameterizedCoefficient)(parameters, time) =
    coefficient.evaluate(parameters, time)

"""Evaluate a scalar expression to either a number or a function of time."""
function _numerical_scalar(expr::ScalarExpr, parameter_lookup)
    expr.head == :const && return expr.args[1]
    expr.head == :param && return parameter_lookup(expr.args[1])
    expr.head == :neg && return _time_unary(-, _numerical_scalar(expr.args[1], parameter_lookup))
    expr.head == :conj && return _time_unary(conj, _numerical_scalar(expr.args[1], parameter_lookup))
    expr.head == :cos && return _time_unary(cos, _numerical_scalar(expr.args[1], parameter_lookup))
    expr.head == :sin && return _time_unary(sin, _numerical_scalar(expr.args[1], parameter_lookup))
    if expr.head in (:+, :*)
        operation = expr.head == :+ ? (+) : (*)
        value = _numerical_scalar(expr.args[1], parameter_lookup)
        for arg in expr.args[2:end]
            value = _time_binary(operation, value, _numerical_scalar(arg, parameter_lookup))
        end
        return value
    elseif expr.head == :/
        return _time_binary(/,
            _numerical_scalar(expr.args[1], parameter_lookup),
            _numerical_scalar(expr.args[2], parameter_lookup))
    elseif expr.head == :^
        return _time_unary(value -> value ^ expr.args[2],
            _numerical_scalar(expr.args[1], parameter_lookup))
    end
    error("Unknown scalar expression head $(expr.head).")
end

_time_unary(operation, value::_ParameterizedCoefficient) =
    _ParameterizedCoefficient((parameters, time) ->
        operation(value(parameters, time)))
_time_unary(operation, value::Function) = time -> operation(value(time))
_time_unary(operation, value) = operation(value)

function _time_binary(operation, left, right)
    if left isa _ParameterizedCoefficient || right isa _ParameterizedCoefficient
        return _ParameterizedCoefficient((parameters, time) ->
            operation(
                _coefficient_at(left, parameters, time),
                _coefficient_at(right, parameters, time),
            ))
    end
    left isa Function || right isa Function || return operation(left, right)
    time -> operation(_at(left, time), _at(right, time))
end

_coefficient_at(value::_ParameterizedCoefficient, parameters, time) =
    value(parameters, time)
_coefficient_at(value::Function, parameters, time) = value(time)
_coefficient_at(value, parameters, time) = value

"""Compile an operator expression using fixed operators and scalar coefficients."""
function _numerical_expression(expr::OperatorExpr, operator_lookup, parameter_lookup)
    terms = _operator_terms(expr, operator_lookup, parameter_lookup)
    static = nothing
    dynamic = Tuple{Any,Function}[]
    for (coefficient, operator) in terms
        if coefficient isa Union{Function,_ParameterizedCoefficient}
            operator isa QuantumObject ||
                error("Time-dependent expressions require QuantumObject operators.")
            push!(dynamic, (operator,
                (parameters, time) ->
                    _coefficient_at(coefficient, parameters, time)))
        else
            term = coefficient * operator
            static = static === nothing ? term : static + term
        end
    end
    isempty(dynamic) && return static
    evolution = QobjEvo(Tuple(dynamic))
    static === nothing ? evolution : static + evolution
end

function _operator_terms(expr::OperatorExpr, operator_lookup, parameter_lookup)
    if expr.head == :op
        return [(1, operator_lookup(expr.args[1].path))]
    elseif expr.head == :+
        return reduce(vcat,
            (_operator_terms(arg, operator_lookup, parameter_lookup) for arg in expr.args))
    elseif expr.head == :*
        terms = [(1, nothing)]
        for arg in expr.args
            if arg isa ScalarExpr
                value = _numerical_scalar(arg, parameter_lookup)
                terms = [(_time_binary(*, coefficient, value), operator)
                    for (coefficient, operator) in terms]
            elseif arg isa Number
                terms = [(coefficient * arg, operator) for (coefficient, operator) in terms]
            else
                terms = _multiply_terms(
                    terms,
                    _operator_terms(arg, operator_lookup, parameter_lookup),
                )
            end
        end
        return terms
    elseif expr.head == :^
        exponent = expr.args[2]
        exponent isa Integer && exponent >= 1 ||
            error("Operator expression powers must be positive integers.")
        base = _operator_terms(expr.args[1], operator_lookup, parameter_lookup)
        terms = base
        for _ in 2:exponent
            terms = _multiply_terms(terms, base)
        end
        return terms
    elseif expr.head == :adjoint
        return [(_time_unary(conj, coefficient), adjoint(operator))
            for (coefficient, operator) in
                _operator_terms(expr.args[1], operator_lookup, parameter_lookup)]
    end
    error("Unknown operator expression head $(expr.head).")
end

function _multiply_terms(left_terms, right_terms)
    [(
        _time_binary(*, left_coefficient, right_coefficient),
        left_operator === nothing ? right_operator :
            right_operator === nothing ? left_operator : left_operator * right_operator,
    ) for (left_coefficient, left_operator) in left_terms
      for (right_coefficient, right_operator) in right_terms]
end

_param_parts(path::Tuple{Vararg{Symbol}}) = path
_param_parts(path::Symbol) = (path,)
_param_parts(path) = hasproperty(path, :path) ? _param_parts(getproperty(path, :path)) : path.parts

"""Return a catalog-relative path, accepting canonical and expression paths."""
function _canonical_reference(path, catalog::Symbol)
    parts = _param_parts(path)
    if length(parts) >= 4 && parts[1] == :components && parts[3] == catalog
        return (parts[2], parts[4:end]...)
    elseif length(parts) >= 2 && parts[1] == catalog
        return parts[2:end]
    end
    return parts
end

if isdefined(@__MODULE__, :DeviceParameter)
    import Base: *

    *(parameter::DeviceParameter, expr::OperatorExpr) = param(parameter) * expr
    *(expr::OperatorExpr, parameter::DeviceParameter) = expr * param(parameter)
end
