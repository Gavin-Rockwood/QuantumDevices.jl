"""A named upstream recipe, its physical owners, and the consuming recipe scope."""
struct InputBinding{O <: AbstractVector{Symbol}, N <: AbstractDict}
    name::Symbol
    recipe::Union{Recipe, Symbol}
    owners::O
    scope::Symbol # identifies symbolic operators bound to this consuming node
    operator_names::N
end

struct InputOperators
    binding::InputBinding
    owner::Union{Nothing, Symbol} # nothing until a grouped input selects a component
end

"""Internal implementation of the public `Recipe(device, components; kwargs...)` constructor."""
struct CompositeRecipe{
    C <: AbstractVector{Symbol},
    S <: AbstractDict,
    P <: AbstractDict,
    I <: AbstractVector,
} <: Recipe
    components::C
    inputs::NamedTuple
    states::S # stable zero-based index or occupation-label requests keyed by alias
    parameters::P
    interactions::I
    max_dimension::Int
    scope::Symbol # local operator binding identity, refreshed when copying
    dressing::Bool
    dressing_options::NamedTuple

    # Direct constructors validate options; contextual capture validates mutable declarations.
    function CompositeRecipe(components::C, inputs::NamedTuple, states::S, parameters::P,
                             interactions::I, max_dimension::Int, scope::Symbol,
                             dressing, dressing_options) where {C <: AbstractVector{Symbol}, S <: AbstractDict, P <: AbstractDict, I <: AbstractVector}
        dressing_options = validate_recipe_dressing(dressing, dressing_options)
        new{C, S, P, I}(components, inputs, states, parameters, interactions, max_dimension, scope, dressing, dressing_options)
    end
end

# Preserve the legacy positional constructor defaults.
CompositeRecipe(components, inputs, states, parameters, interactions, max_dimension, scope) =
    CompositeRecipe(
        components,
        inputs,
        states,
        parameters,
        interactions,
        max_dimension,
        scope,
        true,
        (;),
    )
