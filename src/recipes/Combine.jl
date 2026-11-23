"""A named leaf-recipe input and its ordered, zero-based retained eigenstate indices."""
struct Input{S <: AbstractVector{Int}}
    name::Symbol
    recipe::Union{Recipe, Symbol}
    states::S # zero-based eigenstate indices
end

"""Compatibility recipe composing disjoint named inputs with registered interactions."""
Base.@kwdef struct Combine{I <: AbstractVector, X <: AbstractVector, P <: AbstractDict, M <: AbstractDict} <: Recipe
    inputs::I
    interactions::X = Symbol[]
    parameters::P = Dict{Num, Number}()
    max_dimension::Int = 1000
    dressing::Bool = true
    dressing_options::NamedTuple = (;)
    control_maps::M = Dict{Symbol,AbstractControlMap}()

    # Direct constructors preserve options; capture validates the complete input structure.
    function Combine(inputs::I, interactions::X, parameters::P, max_dimension::Int,
                     dressing, dressing_options, control_maps::M) where {I <: AbstractVector, X <: AbstractVector, P <: AbstractDict, M <: AbstractDict}
        dressing_options = validate_recipe_dressing(dressing, dressing_options)
        new{I, X, P, M}(inputs, interactions, parameters, max_dimension, dressing, dressing_options, control_maps)
    end
end

# Preserve the legacy positional constructor defaults.
Combine(inputs, interactions, parameters, max_dimension) =
    Combine(inputs, interactions, parameters, max_dimension, true, (;))

Combine(inputs,interactions,parameters,max_dimension,dressing,dressing_options) = Combine(inputs,interactions,parameters,max_dimension,dressing,dressing_options,Dict{Symbol,AbstractControlMap}())
