"""A symbolic coupling. Supply numerical inputs through recipe/build parameters."""
struct Interaction{V <: AbstractDict, M <: NamedTuple}
    name::Symbol
    expression::Sym
    variables::V
    metadata::M
end

function Interaction(name, expression; metadata=(;), parameters=nothing)
    isnothing(parameters) || throw(ArgumentError("interaction defaults were removed; put values in recipe.overrides or build(...; parameters=...)"))
    expression = asym(expression)
    variables = Dict{String,Num}(string(v)=>v for v in refs(expression,:param))
    all(is_parameter,values(variables)) || throw(ArgumentError("interaction inputs must be atomic Symbolics variables"))
    return Interaction(Symbol(name),expression,variables,metadata)
end
Interaction(name,expression,defaults;kwargs...) = throw(ArgumentError("interaction defaults were removed; put values in recipe.overrides or build(...; parameters=...)"))

function validate_interaction_inputs(interaction)
    expected = Dict{String,Num}(string(v)=>v for v in refs(interaction.expression,:param))
    isequal(expected,interaction.variables) || throw(ArgumentError("interaction $(interaction.name) variables do not match its expression"))
    return nothing
end
