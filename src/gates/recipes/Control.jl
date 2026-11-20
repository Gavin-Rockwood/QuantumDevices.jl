"""Bind a named pulse's full value to a Hamiltonian parameter (including its offset)."""
struct Control{P <: AbstractPulse}
    name::Symbol
    target::Num
    pulse::P
    function Control(name::Symbol, target, pulse::P) where {P <: AbstractPulse}
        is_parameter(target) || throw(ArgumentError("control target must be an atomic Symbolics variable"))
        return new{P}(name, target, pulse)
    end
end

