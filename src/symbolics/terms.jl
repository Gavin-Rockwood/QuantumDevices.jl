# Each term is a scalar expression and an ordered (reference, adjointed) word.
function terms(symbolic::Sym)
    expression = symbolic.expr
    kind = expression[1]
    isempty(refs(symbolic, :op)) && return [(Any[], symbolic)]
    kind == :op && return [(Any[(expression[2], false)], val(1))]
    if kind == :adjoint
        return [
            (Any[(ref, !adj) for (ref, adj) in reverse(word)], adjoint(coefficient)) for
            (word, coefficient) in terms(expression[2])
        ]
    elseif kind == :power
        return foldl(
            multerms,
            (terms(expression[2]) for _ = 1:expression[3]);
            init = [(Any[], val(1))],
        )
    end
    operation = expression[2]
    left, right = expression[3:4]
    operation === (+) && return vcat(terms(left), terms(right))
    operation === (-) && return vcat(terms(left), [(word, -c) for (word, c) in terms(right)])
    operation === (*) && return multerms(terms(left), terms(right))
    if operation === (/) && isempty(refs(right, :op))
        return [(word, coefficient / right) for (word, coefficient) in terms(left)]
    end
    return throw(ArgumentError("unsupported operator expression"))
end

function multerms(left, right)
    return [
        (vcat(left_word, right_word), left_coefficient * right_coefficient) for
        (left_word, left_coefficient) in left for (right_word, right_coefficient) in right
    ]
end

physical_ref(ref::OperatorRef) = ref
physical_ref(ref::BoundOperatorRef) = ref.physical
function physical_expression(expression::Sym)
    expression.expr[1] == :op && return Sym((:op, physical_ref(expression.expr[2])))
    return Sym(Tuple(x isa Sym ? physical_expression(x) : x for x in expression.expr))
end

