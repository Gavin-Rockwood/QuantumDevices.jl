const _DISPLAY_PREC_SUM = 10
const _DISPLAY_PREC_PRODUCT = 20
const _DISPLAY_PREC_PREFIX = 30
const _DISPLAY_PREC_POWER = 40
const _DISPLAY_PREC_ATOM = 50

function _expression_text(io::IO, expr::ScalarExpr, parent_precedence::Int = 0)
    return _expression_text_impl(io, expr, parent_precedence)
end

function _expression_text(io::IO, expr::OperatorExpr, parent_precedence::Int = 0)
    return _expression_text_impl(io, expr, parent_precedence)
end

function _expression_text_impl(io::IO, expr::Expr, parent_precedence::Int)
    head = expr.head
    if head == :const
        return sprint(show, expr.args[1]; context = IOContext(io, :compact => true))
    elseif head == :param
        parameter = expr.args[1]
        path = parameter isa DeviceParameter ? parameter.path : parameter
        return _display_path(path)
    elseif head == :op
        return _display_operator_path(expr.args[1].path)
    elseif head == :+
        text = _display_sum_text(io, expr.args)
        return _display_parenthesize(text, _DISPLAY_PREC_SUM, parent_precedence)
    elseif head == :*
        text = _display_product_text(io, expr.args)
        return _display_parenthesize(text, _DISPLAY_PREC_PRODUCT, parent_precedence)
    elseif head == :/
        numerator = _expression_text(io, expr.args[1], _DISPLAY_PREC_PRODUCT)
        denominator = _expression_text(io, expr.args[2], _DISPLAY_PREC_PRODUCT + 1)
        text = string(numerator, " / ", denominator)
        return _display_parenthesize(text, _DISPLAY_PREC_PRODUCT, parent_precedence)
    elseif head == :neg
        text = string("-", _expression_text(io, expr.args[1], _DISPLAY_PREC_PREFIX))
        return _display_parenthesize(text, _DISPLAY_PREC_PREFIX, parent_precedence)
    elseif head == :cos || head == :sin
        text = string(head, "(", _expression_text(io, expr.args[1]), ")")
        return _display_parenthesize(text, _DISPLAY_PREC_PREFIX, parent_precedence)
    elseif head == :^
        base = _expression_text(io, expr.args[1], _DISPLAY_PREC_POWER)
        exponent = expr.args[2]
        suffix = _display_exponent(io, exponent)
        text = string(base, suffix)
        return _display_parenthesize(text, _DISPLAY_PREC_POWER, parent_precedence)
    elseif head == :adjoint
        operand = _expression_text(io, expr.args[1], _DISPLAY_PREC_POWER)
        text = _display_unicode(io) ? string(operand, "†") : string("adjoint(", operand, ")")
        return _display_parenthesize(text, _DISPLAY_PREC_POWER, parent_precedence)
    end
    return string(head, "(", join(string.(expr.args), ", "), ")")
end

function _display_sum_text(io::IO, args)
    isempty(args) && return "0"
    text = _expression_text(io, args[1], _DISPLAY_PREC_SUM)
    for arg in args[2:end]
        negative = _display_negative_term(arg)
        if negative === nothing
            text *= " + " * _expression_text(io, arg, _DISPLAY_PREC_SUM)
        else
            text *= " - " * _expression_text(io, negative, _DISPLAY_PREC_SUM + 1)
        end
    end
    return text
end

function _display_negative_term(expr)
    if expr isa ScalarExpr && expr.head == :neg
        return expr.args[1]
    elseif expr isa Expr && expr.head == :* && !isempty(expr.args)
        first_factor = expr.args[1]
        if first_factor isa ScalarExpr && first_factor.head == :const &&
                first_factor.args[1] == -1
            remaining = expr.args[2:end]
            isempty(remaining) && return constant(1)
            length(remaining) == 1 && return remaining[1]
            return expr isa ScalarExpr ?
                ScalarExpr(:*, Tuple(remaining)) :
                OperatorExpr(:*, Tuple(remaining))
        end
    end
    return nothing
end

function _display_product_text(io::IO, args)
    factors = collect(args)
    if !isempty(factors)
        final = factors[end]
        if final isa ScalarExpr && final.head == :/ &&
                _display_is_one(final.args[1])
            prefix = factors[1:end-1]
            product = join(
                (_expression_text(io, factor, _DISPLAY_PREC_PRODUCT) for factor in prefix),
                _display_symbol(io, " × ", " * "),
            )
            denominator = _expression_text(io, final.args[2], _DISPLAY_PREC_PRODUCT + 1)
            return isempty(product) ? string("1 / ", denominator) : string(product, " / ", denominator)
        end
    end
    return join(
        (_expression_text(io, factor, _DISPLAY_PREC_PRODUCT) for factor in factors),
        _display_symbol(io, " × ", " * "),
    )
end

_display_is_one(expr::ScalarExpr) =
    expr.head == :const && length(expr.args) == 1 && expr.args[1] == 1

function _display_exponent(io::IO, exponent)
    if _display_unicode(io) && exponent isa Integer
        superscripts = Dict(
            '-' => '⁻', '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³',
            '4' => '⁴', '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
        )
        return join(get(superscripts, character, character) for character in string(exponent))
    end
    return string("^", exponent)
end

_display_parenthesize(text, precedence, parent_precedence) =
    precedence < parent_precedence ? string("(", text, ")") : text

function Base.show(io::IO, ref::OperatorRef)
    print(io, "OperatorRef(", _display_operator_path(ref.path), ")")
end

Base.show(io::IO, expr::Expr) = print(io, _expression_text(io, expr))
Base.show(io::IO, ::MIME"text/plain", expr::Expr) = print(io, _expression_text(io, expr))
