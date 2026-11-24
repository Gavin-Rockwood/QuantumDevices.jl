"""A symbolic coupling. Declared defaults are captured by recipes; recipe/build parameters override them."""
struct Interaction{V <: AbstractDict, D <: AbstractDict, M <: NamedTuple}
    name::Symbol
    expression::Sym
    variables::V
    defaults::D
    metadata::M
end

function Interaction(name, expression; metadata=(;), parameters=Dict(), defaults=parameters)
    expression = asym(expression)
    variables = Dict{String,Num}(string(v)=>v for v in refs(expression,:param))
    all(is_parameter,values(variables)) || throw(ArgumentError("interaction inputs must be atomic Symbolics variables"))
    all(v -> any(isequal(v),values(variables)),keys(defaults)) || throw(ArgumentError("interaction defaults refer to unknown variables"))
    declared = Dict{Num,Union{Missing,Number}}(v=>get(defaults,v,missing) for v in values(variables))
    all(x->ismissing(x) || concrete_parameter_value(x),values(declared)) ||
        throw(ArgumentError("interaction defaults must be finite concrete numbers or missing"))
    return Interaction(Symbol(name),expression,variables,declared,metadata)
end
Interaction(name,expression,defaults;kwargs...) = Interaction(name,expression;defaults,kwargs...)
# Compatibility for internal declarations without numerical defaults.
Interaction(name,expression,variables,metadata) = Interaction(name,expression,variables,
    Dict{Num,Union{Missing,Number}}(v=>missing for v in values(variables)),metadata)

function validate_interaction_inputs(interaction)
    expected = Dict{String,Num}(string(v)=>v for v in refs(interaction.expression,:param))
    isequal(expected,interaction.variables) || throw(ArgumentError("interaction $(interaction.name) variables do not match its expression"))
    Set(keys(interaction.defaults)) == Set(values(interaction.variables)) ||
        throw(ArgumentError("interaction defaults do not match its declared variables"))
    all(x->ismissing(x) || concrete_parameter_value(x),values(interaction.defaults)) ||
        throw(ArgumentError("interaction defaults must be finite concrete numbers or missing"))
    return nothing
end
