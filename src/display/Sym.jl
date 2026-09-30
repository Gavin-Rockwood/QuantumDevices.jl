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
