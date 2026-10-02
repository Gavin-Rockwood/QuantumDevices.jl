# ============================================================
# Symbolic expression
# ============================================================

"""
    Sym(expr)

Symbolic operator/scalar expression used to construct Hamiltonians. Prefer
[`op`](@ref), [`param`](@ref), [`val`](@ref), and [`call`](@ref) to constructing
its internal tuple representation directly. Arithmetic, integer powers, `sin`,
`cos`, `exp`, and `sqrt` compose expressions without evaluating them.
Operator multiplication preserves order; names are resolved by [`numerical`](@ref).
"""
struct Sym
    expr
end

"""
    op(name::Symbol)

Return a symbolic operator reference. Component-local names are such as `:charge`;
model names are prefixed, such as `:t_charge`. Referenced names must exist in the
operator container supplied to [`numerical`](@ref).

```jldoctest
julia> op(:x) isa Sym
true
```
"""
op(name::Symbol) = Sym((:op, name))
"""
    param(name::Symbol)

Return a symbolic scalar parameter reference. Model component parameters use
`componentname_localname`; coupling and gate parameters use their supplied names.
Missing parameters in low-level [`numerical`](@ref) remain dynamic; model and gate
construction reject unresolved required values.
"""
param(name::Symbol) = Sym((:param, name))
"""
    val(x::Number)

Wrap a numeric constant in a [`Sym`](@ref). Constants are independent of parameter
containers and time. Useful for a zero interaction: `val(0)`.
"""
val(x::Number) = Sym((:val, x))
_symbolic_argument(x) = x
_symbolic_argument(x::Number) = val(x)
"""
    call(f, args...)

Construct a symbolic function application. Numeric arguments are wrapped as
[`val`](@ref). Operator-valued expressions must fit the supported sum/product
compiler; general functions are intended for scalar coefficients. Persistence
supports only its explicit built-in operation allowlist, not arbitrary functions.
"""
call(f, args...) = Sym((:call, f, map(_symbolic_argument, args)...))
