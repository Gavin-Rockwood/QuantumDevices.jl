"""
    Envelope(f, parameters=(;))

Custom dimensionless envelope with explicit named parameters. `f(p, t, duration)`
receives local pulse time and duration. Parameters are exposed recursively to
`getpath`, `setpath`, and calibration. Callable persistence uses JLD2 and requires
its definition to be available when loading; built-in shapes have no such requirement.
"""
struct Envelope{F,P<:NamedTuple} <: AbstractEnvelope
    callable::F
    parameters::P
    function Envelope(f::F, p::P) where {F,P<:NamedTuple}
        applicable(f, p, 0.0, 1.0) || throw(ArgumentError("Envelope must accept f(p, t, duration)"))
        new{F,P}(f, p)
    end
end
Envelope(f) = Envelope(f, (;))
