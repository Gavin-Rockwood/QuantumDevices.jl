# One expression renderer supplies both compact text and notebook mathematics.
const DISPLAY_GREEK = Dict(
    "ν" => "nu",
    "ϕ" => "varphi",
    "φ" => "varphi",
    "α" => "alpha",
    "β" => "beta",
    "γ" => "gamma",
    "δ" => "delta",
    "Δ" => "Delta",
    "ω" => "omega",
    "Ω" => "Omega",
    "λ" => "lambda",
    "π" => "pi",
    "θ" => "theta",
    "μ" => "mu",
    "χ" => "chi",
)
function texescape(value)
    substitutions = Dict(
        '\\' => "\\textbackslash{}",
        '_' => "\\_",
        '{' => "\\{",
        '}' => "\\}",
        '%' => "\\%",
        '&' => "\\&",
        '#' => "\\#",
        '$' => raw"\$",
        '^' => "\\textasciicircum{}",
        '~' => "\\textasciitilde{}",
    )
    return join(get(substitutions, c, string(c)) for c in string(value))
end
identifier(value) =
    haskey(DISPLAY_GREEK, string(value)) ? "\\" * DISPLAY_GREEK[string(value)] :
    "\\mathrm{" * texescape(value) * "}"
function parameter_label(expression::Number, labels=Dict())
    isempty(labels) && return string(Latexify.latexify(expression; env=:raw))
    replacements = Dict(Symbolics.toexpr(variable)=>latexstring(label_value(label))
        for (variable,label) in labels if variable isa Num)
    rendered = scalar_display_labels(Symbolics.toexpr(expression), replacements)
    return string(Latexify.latexify(rendered; env=:raw))
end
scalar_display_labels(expression, labels) = get(labels, expression, expression)
function scalar_display_labels(expression::Expr, labels)
    haskey(labels, expression) && return labels[expression]
    return Expr(expression.head, (scalar_display_labels(argument,labels) for argument in expression.args)...)
end

function operator_label(ref::OperatorRef)
    symbol = get(
        Dict(
            :charge => "n",
            :a => "a",
            :number => "N",
            :n => "N",
            :jump => "J",
            :x => "\\sigma^{x}",
            :y => "\\sigma^{y}",
            :z => "\\sigma^{z}",
        ),
        ref.name,
        nothing,
    )
    isnothing(symbol) && (symbol = identifier(ref.name))
    return "\\hat{" * symbol * "}_{" * identifier(ref.owner) * "}"
end
label_value(value::LaTeXString) = replace(String(value), r"^\$|\$$" => "")
label_value(value::AbstractString) = "\\mathrm{" * texescape(value) * "}"
label_value(value) = throw(ArgumentError("display labels must be String or LaTeXString values"))
function mathlabel(object, labels)
    haskey(labels, object) && return label_value(labels[object])
    object isa Num && return parameter_label(object)
    object isa OperatorRef && return operator_label(object)
    if object isa BoundOperatorRef
        annotation =
            object.input == object.physical.owner ? "" : "^{[" * identifier(object.input) * "]}"
        return mathlabel(object.physical, labels) * annotation
    end
    return throw(ArgumentError("unsupported mathematical label $(typeof(object))"))
end

function number_math(value::Number, tex)
    if value isa Complex
        realpart, imaginary = real(value), imag(value)
        unit = tex ? "\\mathrm{i}" : "im"
        return "(" *
               number_math(realpart, tex) *
               (imaginary < 0 ? " - " : " + ") *
               number_math(abs(imaginary), tex) *
               (tex ? "\\," : "*") *
               unit *
               ")"
    elseif value isa Rational
        return tex ?
               "\\frac{" * string(numerator(value)) * "}{" * string(denominator(value)) * "}" :
               "(" * string(numerator(value)) * "/" * string(denominator(value)) * ")"
    end
    return string(value)
end
paren(s, tex) = tex ? "\\left(" * s * "\\right)" : "(" * s * ")"
# Higher precedence binds more tightly; equal-precedence right operands need grouping
# for subtraction/division and for preserving a nested right-hand arithmetic tree.
function render_expression(s::Sym, tex::Bool, labels = Dict(); parent = 0)
    expression = s.expr
    kind = expression[1]
    precedence = 100
    if kind == :val
        result = (expression[2] isa Num || expression[2] isa Complex{<:Num}) && tex ? (haskey(labels, expression[2]) ? mathlabel(expression[2],labels) : parameter_label(expression[2],labels)) : number_math(expression[2], tex)
        !(expression[2] isa Num) && expression[2] isa Real && expression[2] < 0 && (precedence = 25)
    elseif kind == :op
        result = tex ? mathlabel(expression[2], labels) : sprint(show, expression[2])
    elseif kind == :power
        precedence = 30
        base = render_expression(expression[2], tex, labels; parent = 31)
        result = base * (tex ? "^{" * string(expression[3]) * "}" : "^" * string(expression[3]))
    elseif kind == :adjoint
        precedence = 40
        operand = render_expression(expression[2], tex, labels; parent = 41)
        result = if isempty(refs(expression[2], :op))
            tex ? "\\overline{" * render_expression(expression[2], tex, labels) * "}" :
            "conj(" * render_expression(expression[2], tex, labels) * ")"
        else
            operand * (tex ? "^{\\dagger}" : "′")
        end
    elseif kind == :call
        operation = expression[2]
        args = expression[3:end]
        if operation in (+, -, *, /) && length(args) == 2
            precedence = operation in (+, -) ? 10 : 20
            left = render_expression(args[1], tex, labels; parent = precedence)
            right = render_expression(args[2], tex, labels; parent = precedence + 1)
            if operation === (/) && tex
                result =
                    "\\frac{" *
                    render_expression(args[1], tex, labels) *
                    "}{" *
                    render_expression(args[2], tex, labels) *
                    "}"
            elseif operation === (*) && isequal(args[1].expr, (:val, -1))
                result = "-" * render_expression(args[2], tex, labels; parent = 21)
            else
                separator =
                    operation === (*) ? (tex ? " \\, " : " * ") :
                    " " * string(nameof(operation)) * " "
                result = left * separator * right
            end
        elseif operation === sqrt && length(args) == 1 && tex
            result = "\\sqrt{" * render_expression(only(args), tex, labels) * "}"
        elseif operation in (cospi, sinpi) && length(args) == 1 && tex
            result =
                (operation === cospi ? "\\cos" : "\\sin") *
                paren("\\pi \\, " * render_expression(only(args), tex, labels; parent = 21), tex)
        else
            name =
                tex ? "\\operatorname{" * texescape(nameof(operation)) * "}" :
                string(nameof(operation))
            result =
                name * paren(join((render_expression(a, tex, labels) for a in args), ", "), tex)
        end
    else
        throw(ArgumentError("unsupported symbolic display node $kind"))
    end
    return precedence < parent ? paren(result, tex) : result
end

const MathDisplayObject = Union{Sym, OperatorRef, BoundOperatorRef}
function mathbody(object::MathDisplayObject, labels)
    object isa Sym && return render_expression(object, true, labels)
    return mathlabel(object, labels)
end
