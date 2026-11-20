Base.show(io::IO, ref::OperatorRef) = print(io, ref.owner, ".", ref.name)

for f in (:+, :-, :*, :/)
    @eval begin
        Base.$f(a::Sym, b::Sym) = isempty(refs(a, :op)) && isempty(refs(b, :op)) ? val($f(scalar_expression(a), scalar_expression(b))) : Sym((:call, $f, a, b))
        Base.$f(a::Sym, b::Number) = $f(a, asym(b))
        Base.$f(a::Number, b::Sym) = $f(asym(a), b)
    end
end

Base.:-(a::Sym) = val(-1) * a
Base.:^(a::Sym, n::Integer) = n >= 0 ? Sym((:power, a, n)) : throw(ArgumentError("negative power"))
Base.literal_pow(::typeof(^), a::Sym, ::Val{N}) where {N} = a ^ N
Base.adjoint(a::Sym) = isempty(refs(a, :op)) ? val(scalar_conjugate(scalar_expression(a))) : Sym((:adjoint, a))

# These functions are scalar coefficients, never operator matrix functions.
for f in (:cospi, :sinpi, :sqrt)
    @eval function Base.$f(a::Sym)
        isempty(refs(a, :op)) || throw(ArgumentError("scalar functions cannot contain operators"))
        return val($f(scalar_expression(a)))
    end
end

"""
    latex(object; labels=Dict())

Return notebook mathematics without altering symbolic identities. Label keys are
Symbolics variables, OperatorRef, or bound operator identities; LaTeXString values are raw
math, while ordinary strings are escaped labels. Bound operators also accept a
physical OperatorRef fallback. Input annotations are omitted when the alias equals
the physical owner.
"""
function latex(object::MathDisplayObject; labels = Dict())
    return latex_object(object, labels)
end

function latex_object(object, labels)
    all(k -> k isa Union{Num, OperatorRef, BoundOperatorRef}, keys(labels)) ||
        throw(ArgumentError("label keys must be parameter or operator identities"))
    return latexstring(mathbody(object, labels))
end

Base.show(io::IO, ::MIME"text/latex", object::MathDisplayObject) = print(io, latex(object))
Base.show(io::IO, s::Sym) = print(io, render_expression(s, false))

Base.show(io::IO, ref::BoundOperatorRef) = print(io, ref.input, ".operators.", ref.physical)

# A package helper only: Symbolics owns Num printing globally.
latex(object::Num; labels=Dict()) = latex_object(val(object), labels)

# Keep explicit conjugation even for real-declared variables: runtime controls may
# supply complex drive amplitudes. Complex scalar leaves need the whole conjugate.
scalar_conjugate(x::Complex) = scalar_conjugate(real(x)) - im*scalar_conjugate(imag(x))
scalar_conjugate(x::Num) = Num(Symbolics.term(conj, x))
scalar_conjugate(x::Number) = conj(x)
