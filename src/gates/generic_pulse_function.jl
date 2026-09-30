"""
    GenericPulseFunction(callable, parameters::NamedTuple)

Construct a pulse with `callable(p, t)` and named parameters. Evaluation supplies
`p.duration`; do not store that field yourself. A callable must accept the two
arguments and return a finite scalar number. The parameters participate in
[`parameters`](@ref) discovery and [`setpath`](@ref) updates.

Persistence uses JLD2 for the callable. Put portable callable definitions in a
module loaded before restoring a bundle; arbitrary interactive closures are not
portable source archives.
"""
struct GenericPulseFunction{F,P<:NamedTuple} <: AbstractParameterizedPulse
    callable::F
    parameters::P

    function GenericPulseFunction(callable::F, parameters::P) where {F,P<:NamedTuple}
        _check_pulse_parameters(parameters)
        applicable(callable, merge(parameters, (; duration = 1.0)), 0.0) ||
            throw(ArgumentError("Generic pulse must be callable as f(p, t)"))
        new{F,P}(callable, parameters)
    end
end

pulse_function(pulse::GenericPulseFunction) = pulse.callable
