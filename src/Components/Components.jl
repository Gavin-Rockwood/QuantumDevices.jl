"""
    AbstractComponentSpec

Base type for component specifications that provide parameters, operator
catalogs, and optionally a symbolic Hamiltonian.
"""
abstract type AbstractComponentSpec end

function _component_hamiltonian(hamiltonian, extra_terms)
    hamiltonian === nothing && return extra_terms
    extra_terms === nothing && return hamiltonian
    hamiltonian + extra_terms
end

"""
    numerical_operator(spec, operator; kwargs...)

Evaluate one operator from a component specification's operator catalog.
"""
function numerical_operator(spec::AbstractComponentSpec, operator, args...; kwargs...)
    if !(operator in keys(spec.operators))
        error("Operator $operator not in Spec Operators")
    end
    return spec.operators[operator](args...; kwargs...)
end

include("component.jl")
include("repository/qubit.jl")
include("repository/generic.jl")
include("repository/transmon.jl")
include("repository/flux_tunable_transmon.jl")
include("repository/resonator.jl")
