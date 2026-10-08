const Container = Union{NamedTuple, AbstractDict}

parameter_keys(x) = Set{Symbol}()
function parameter_keys(x::Sym)
    kind, args... = x.expr
    kind === :param && return Set([only(args)])
    kind === :op && return Set{Symbol}()
    kind === :val && return Set{Symbol}()
    kind === :call || throw(ArgumentError("Unknown symbolic expression: $(x.expr)"))
    return union(Set{Symbol}(), (parameter_keys(a) for a in args[2:end])...)
end

_has_operator(x) = false
function _has_operator(x::Sym)
    x.expr[1] === :op && return true
    x.expr[1] === :param && return false
    x.expr[1] === :val && return false
    x.expr[1] === :call || throw(ArgumentError("Unknown symbolic expression: $(x.expr)"))
    return any(_has_operator, x.expr[3:end])
end

numerical_static(x, params, operators) = x
function numerical_static(x::Sym, params::Container, operators::Container)
    e = x.expr
    e[1] === :param && return _get(params, e[2])
    e[1] === :op && return _get(operators, e[2])
    e[1] === :val && return e[2]
    e[1] === :call || throw(ArgumentError("Unknown symbolic expression: $e"))
    return e[2](map(a -> numerical_static(a, params, operators), e[3:end])...)
end

_resolve_param(x, t) = x
_resolve_param(f::Function, t) = f(t)

function scalar_function(x, fixed_params::Container)
    x isa Sym || return (p, t) -> x
    e = x.expr
    if e[1] === :param
        key = e[2]
        if _has(fixed_params, key)
            value = _get(fixed_params, key)
            return (p, t) -> _resolve_param(value, t)
        end
        return (p, t) -> _resolve_param(_get(p, key), t)
    elseif e[1] === :val
        return (p, t) -> e[2]
    elseif e[1] === :call
        fs = map(a -> scalar_function(a, fixed_params), e[3:end])
        return (p, t) -> e[2]((f(p, t) for f in fs)...)
    end
    throw(ArgumentError("Expected scalar expression, got $x"))
end

# Keep ordered operator words intact until the caller chooses a basis/projection.
_multiply_terms(A, B) = [(vcat(wa, wb), ca * cb) for (wa, ca) in A for (wb, cb) in B]
function _symbolic_terms(x)
    !_has_operator(x) && return [(Symbol[], x)]
    e = x.expr
    e[1] === :op && return [([e[2]], 1)]
    f, args = e[2], e[3:end]
    if f === (+)
        return reduce(vcat, map(_symbolic_terms, args))
    elseif f === (-)
        A = _symbolic_terms(first(args))
        length(args) == 1 && return [(w, -c) for (w, c) in A]
        return vcat(A, [(w, -c) for (w, c) in _symbolic_terms(args[2])])
    elseif f === (*)
        return reduce(_multiply_terms, map(_symbolic_terms, args))
    elseif f === (/)
        _has_operator(args[2]) && throw(ArgumentError("Operator denominators are unsupported"))
        return [(w, c / args[2]) for (w, c) in _symbolic_terms(args[1])]
    elseif f === (^)
        base, n = args
        n = n isa Sym && n.expr[1] === :val ? n.expr[2] : n
        n isa Integer && n >= 0 || throw(ArgumentError("Operator powers must be nonnegative integers"))
        result = [(Symbol[], 1)]
        for _ in 1:n
            result = _multiply_terms(result, _symbolic_terms(base))
        end
        return result
    end
    throw(ArgumentError("Unsupported operator expression: $x"))
end

# Extract multiplicative constants before constructing time-dependent operators.
# Constants inside nonlinear expressions stay in the coefficient callback.
function _constant_coefficient(x, params)
    return all(k -> _has(params, k) && !(_get(params, k) isa Function), parameter_keys(x))
end

function _split_coefficient(x, params)
    _constant_coefficient(x, params) && return scalar_function(x, params)(params, 0.0), nothing
    e = x.expr
    if e[1] === :call && e[2] === (*)
        factor = 1
        residual = nothing
        for arg in e[3:end]
            constant, dynamic = _split_coefficient(arg, params)
            factor *= constant
            dynamic === nothing && continue
            residual = residual === nothing ? dynamic : residual * dynamic
        end
        return factor, residual
    elseif e[1] === :call && e[2] === (/) && _constant_coefficient(e[4], params)
        factor, residual = _split_coefficient(e[3], params)
        return factor / scalar_function(e[4], params)(params, 0.0), residual
    elseif e[1] === :call && e[2] === (-) && length(e) == 3
        factor, residual = _split_coefficient(e[3], params)
        return -factor, residual
    end
    return 1, x
end

function qobjevo_terms(x, operators::Container)
    isempty(operators) && throw(ArgumentError("operators cannot be empty"))
    identity = first(values(operators))^0
    return [(foldl((A, key) -> A * _get(operators, key), word; init = identity), coefficient)
            for (word, coefficient) in _symbolic_terms(x)]
end

# Storage conversion happens at construction, never inside coefficient callbacks.
_numerical_storage(x, ::Bool) = x
_numerical_storage(x::QuantumObject, dense::Bool) =
    dense && !(x.data isa Matrix) ? QuantumObject(Matrix(x.data), x.type, x.dimensions) : x

function build_qobjevo(H, fixed_params::Container, operators::Container; scalar::Number=1, dense::Bool=false)
    result = _evaluate_terms(qobjevo_terms(H, operators), fixed_params; scalar, dense)
    return result isa QobjEvo ? result : QobjEvo(result)
end

function _evaluate_terms(terms, params::Container; fixed=nothing, scalar::Number=1, dense::Bool=false)
    constant = _numerical_storage(fixed, dense) # Cached terms are already scaled.
    dynamic = []
    for (operator, coefficient) in terms
        factor, residual = _split_coefficient(coefficient, params)
        matrix = _numerical_storage((scalar * factor) * operator, dense)
        if residual === nothing
            constant = constant === nothing ? matrix : constant + matrix
        else
            push!(dynamic, (matrix, scalar_function(residual, params)))
        end
    end
    isempty(dynamic) && return constant
    terms_ = constant === nothing ? tuple(dynamic...) : (constant, dynamic...)
    return QobjEvo(terms_)
end
