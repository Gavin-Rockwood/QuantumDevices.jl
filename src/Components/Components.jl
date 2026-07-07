"""
    AbstractComponentSpec

Base type for component specifications that provide parameters, operator
catalogs, and optionally a symbolic Hamiltonian.
"""
abstract type AbstractComponentSpec end

"""
    available_operators(spec)

Return the operator names exposed by a component specification.
"""
available_operators(spec::AbstractComponentSpec) = Set(collect(keys(spec.operators)))

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

include("genarics.jl")
include("component.jl")
include("repository/qubit.jl")
include("repository/generic.jl")
include("repository/transmon.jl")
include("repository/flux_tunable_transmon.jl")
include("repository/resonator.jl")
