"""
    OperatorRef(path)

Symbolic path to an operator in either a component catalog or a local component
catalog.
"""
struct OperatorRef
    path::Tuple{Vararg{Symbol}}
end

abstract type ExprRole end
struct ScalarRole <: ExprRole end
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

const ScalarExpr = Expr{ScalarRole}
const OperatorExpr = Expr{OperatorRole}

constant(value) = ScalarExpr(:const, (value,))

"""
    param(path)

Create a scalar parameter reference for a symbol, tuple path, string path, or
`DeviceParameter`.
"""
param(path) = ScalarExpr(:param, (_expr_path(path),))

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

_expr_path(path::Symbol) = (path,)
_expr_path(path::Tuple{Vararg{Symbol}}) = path
_expr_path(path::AbstractVector{Symbol}) = Tuple(path)
_expr_path(path::AbstractString) = Tuple(Symbol.(split(path, ".")))
_expr_path(path) = hasproperty(path, :path) && hasproperty(path, :default) ?
    path :
    (isdefined(@__MODULE__, :_param_path) ? _param_path(path) : path)

_scalar(expr::ScalarExpr) = expr
_scalar(value::Number) = constant(value)
_scalar(value) = hasproperty(value, :path) && hasproperty(value, :default) ? param(value) : value

import Base: +, -, *, /, ^, adjoint, cos, sin

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
    _participating_components(expr)

Return component IDs explicitly referenced by component-rooted operator paths
inside `expr`.
"""
function _participating_components(expr::OperatorExpr)
    components = Set{Symbol}()
    function collect!(node::OperatorExpr)
        if node.head == :op
            path = node.args[1].path
            length(path) >= 2 && path[1] == :components &&
                push!(components, path[2])
            return
        end
        for arg in node.args
            arg isa OperatorExpr && collect!(arg)
        end
    end
    collect!(expr)
    return components
end

_participating_components(::ScalarExpr) = Set{Symbol}()

"""
    absolute_path(expr, root)

Rewrite relative `param` and `op` references so they are rooted at `root`.
"""
function absolute_path(expr::ScalarExpr, root::Tuple{Vararg{Symbol}})
    if expr.head == :param
        return ScalarExpr(:param, (_absolute_param_path(expr.args[1], root),))
    end
    return ScalarExpr(expr.head, Tuple(_absolute_arg(arg, root) for arg in expr.args))
end

function absolute_path(expr::OperatorExpr, root::Tuple{Vararg{Symbol}})
    if expr.head == :op
        return OperatorExpr(:op, (OperatorRef(_absolute_operator_path(expr.args[1].path, root)),))
    end
    return OperatorExpr(expr.head, Tuple(_absolute_arg(arg, root) for arg in expr.args))
end

absolute_path(expr::Expr, component) =
    absolute_path(expr, (:components, component.id))

_absolute_arg(arg::Expr, root::Tuple{Vararg{Symbol}}) = absolute_path(arg, root)
_absolute_arg(arg, ::Tuple{Vararg{Symbol}}) = arg

function _absolute_operator_path(path::Tuple{Vararg{Symbol}}, root::Tuple{Vararg{Symbol}})
    if length(path) >= 1 && path[1] == :components
        return path
    elseif length(path) >= 1 && path[1] == :operators
        return (root..., path...)
    else
        return (root..., :operators, path...)
    end
end

function _absolute_param_path(path::Tuple{Vararg{Symbol}}, root::Tuple{Vararg{Symbol}})
    if length(path) >= 1 && path[1] == :components
        return path
    elseif length(path) >= 1 && path[1] == :params
        return (root..., path...)
    else
        return (root..., :params, path...)
    end
end

_absolute_param_path(path, root::Tuple{Vararg{Symbol}}) =
    _absolute_param_path(_param_parts(path), root)

_param_parts(path::Tuple{Vararg{Symbol}}) = path
_param_parts(path::Symbol) = (path,)
_param_parts(path) = hasproperty(path, :path) ? _param_parts(getproperty(path, :path)) : path.parts

if isdefined(@__MODULE__, :DeviceParameter)
    import Base: *

    *(parameter::DeviceParameter, expr::OperatorExpr) = param(parameter) * expr
    *(expr::OperatorExpr, parameter::DeviceParameter) = expr * param(parameter)
end
