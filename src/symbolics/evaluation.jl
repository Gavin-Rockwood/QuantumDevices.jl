function children(s::Sym)
    kind = s.expr[1]
    kind == :call && return s.expr[3:end]
    kind in (:power, :adjoint) && return (s.expr[2],)
    return ()
end

function refs(s::Sym, kind)
    s.expr[1] == kind && return Set([s.expr[2]])
    kind == :param && s.expr[1] == :val && return Set(Symbolics.get_variables(s.expr[2]))
    return union(Set{Any}(), (refs(a, kind) for a in children(s))...)
end

"""Evaluate a symbolic tree using parameter and owner-aware operator dictionaries."""
function evaluate(s::Sym, parameters, operators)
    expression = s.expr
    kind = expression[1]
    kind == :val && return scalar_value(expression[2], parameters)
    kind == :op && return operators[expression[2]]
    kind == :adjoint && return adjoint(evaluate(expression[2], parameters, operators))
    kind == :power && return evaluate(expression[2], parameters, operators) ^ expression[3]
    return expression[2]((evaluate(a, parameters, operators) for a in expression[3:end])...)
end

# Scalar arithmetic is delegated entirely to Symbolics; Sym owns operator order.
function scalar_expression(s::Sym)
    isempty(refs(s, :op)) || throw(ArgumentError("expected scalar coefficient"))
    kind = s.expr[1]
    kind == :val && return s.expr[2]
    kind == :power && return scalar_expression(s.expr[2]) ^ s.expr[3]
    kind == :adjoint && return conj(scalar_expression(s.expr[2]))
    return s.expr[2](map(scalar_expression, s.expr[3:end])...)
end

evaluate(s::Number, parameters, operators) = scalar_value(s, parameters)
