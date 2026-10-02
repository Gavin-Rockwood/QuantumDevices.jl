# ============================================================
# Arithmetic
# ============================================================

Base.:+(a::Sym, b::Sym) = call(+, a, b)
Base.:+(a::Sym, b) = call(+, a, b)
Base.:+(a, b::Sym) = call(+, a, b)

Base.:-(a::Sym, b::Sym) = call(-, a, b)
Base.:-(a::Sym, b) = call(-, a, b)
Base.:-(a, b::Sym) = call(-, a, b)
Base.:-(a::Sym) = call(-, a)

Base.:*(a::Sym, b::Sym) = call(*, a, b)
Base.:*(a::Sym, b) = call(*, a, b)
Base.:*(a, b::Sym) = call(*, a, b)

Base.:/(a::Sym, b::Sym) = call(/, a, b)
Base.:/(a::Sym, b) = call(/, a, b)
Base.:/(a, b::Sym) = call(/, a, b)

Base.:^(a::Sym, b) = call(^, a, b)
Base.:^(a::Number, b::Sym) = call(^, a, b)
Base.literal_pow(::typeof(^), a::Sym, ::Val{n}) where {n} = call(^, a, n)

Base.sin(a::Sym) = call(sin, a)
Base.cos(a::Sym) = call(cos, a)
Base.exp(a::Sym) = call(exp, a)
Base.sqrt(a::Sym) = call(sqrt, a)

function Base.show(io::IO, x::Sym)
    e = x.expr
    kind = e[1]
    if kind === :val
        show(io, e[2])
    elseif kind === :op || kind === :param
        print(io, kind, "(")
        show(io, e[2])
        print(io, ")")
    elseif kind === :call
        f, args = e[2], e[3:end]
        if f in (+, -, *, /, ^) && length(args) == 2
            print(io, "(")
            show(io, args[1])
            print(io, " ", f, " ")
            show(io, args[2])
            print(io, ")")
        else
            show(io, f)
            print(io, "(")
            for (i, arg) in enumerate(args)
                i > 1 && print(io, ", ")
                show(io, arg)
            end
            print(io, ")")
        end
    else
        print(io, "Sym(")
        show(io, e)
        print(io, ")")
    end
end

function Base.show(io::IO, ::MIME"text/plain", x::Sym)
    print(io, "Sym: ")
    show(io, x)
end

"""
    numerical(H::Sym, operators[, params])

Evaluate a symbolic expression. Complete scalar parameters produce a static result;
missing or function-valued parameters produce a `QobjEvo`, callable as `(params, t)`.
Function-valued parameters receive time as their single argument.
"""
function numerical(H::Sym, operators::Container, params::Container)
    required = parameter_keys(H)
    dynamic = any(k -> !_has(params, k) || _get(params, k) isa Function, required)
    return dynamic ? build_qobjevo(H, params, operators) : numerical_static(H, params, operators)
end
numerical(H::Sym, operators::Container) = build_qobjevo(H, NamedTuple(), operators)
