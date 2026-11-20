"""An operator identity given by component owner and local operator name."""
struct OperatorRef
    owner::Symbol
    name::Symbol
end

"""An operator supplied by one input in a particular consuming recipe scope."""
struct BoundOperatorRef
    scope::Symbol
    input::Symbol
    physical::OperatorRef
end


"""An owner-aware symbolic arithmetic tree for Hamiltonians."""
struct Sym
    expr::Tuple
end

"""Reference a local operator belonging to `owner`."""
op(owner::Symbol, name::Symbol) = Sym((:op, OperatorRef(owner, name)))
"""Insert a numerical constant into a symbolic expression."""
val(n::Number) = Sym((:val, n))
asym(s::Sym) = s
asym(n::Number) = val(n)
