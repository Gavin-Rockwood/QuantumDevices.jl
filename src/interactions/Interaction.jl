"""A coupling with all used variables and optional explicitly declared defaults.

String declaration keys must be the full symbolic name and map to `(variable, value)`.
Atomic symbolic keys may map directly to values. Variables retain their identities.
"""
struct Interaction{V <: AbstractDict, P <: AbstractDict, M <: NamedTuple}
    name::Symbol
    expression::Sym
    variables::V
    parameters::P
    metadata::M
end

function Interaction(name, expression, defaults=Dict(); parameters=Dict(), metadata=(;))
    expression = asym(expression)
    variables = Dict{String,Num}(string(v)=>v for v in refs(expression, :param))
    all(is_parameter,values(variables)) ||
        throw(ArgumentError("interaction variables must be atomic Symbolics scalar variables"))
    declarations = Dict{String,Number}()
    for (key, declaration) in vcat(collect(pairs(defaults)), collect(pairs(parameters)))
        variable, value = if key isa AbstractString
            declaration isa Tuple && length(declaration) == 2 ||
                throw(ArgumentError("string defaults require (variable, value)"))
            declaration
        else
            key, declaration
        end
        is_parameter(variable) && concrete_parameter_value(value) ||
            throw(ArgumentError("defaults require atomic variables and finite concrete numbers"))
        alias = string(variable)
        key isa AbstractString && key != alias &&
            throw(ArgumentError("default key $key must equal full symbolic name $alias"))
        haskey(variables, alias) && isequal(variables[alias],variable) ||
            throw(ArgumentError("default variable $variable is absent from interaction $name"))
        haskey(declarations,alias) && throw(ArgumentError("duplicate default declaration for $variable"))
        declarations[alias] = value
    end
    return Interaction(Symbol(name), expression, variables, declarations, metadata)
end
function interaction_parameter_values(i::Interaction)
    expected = Dict{String,Num}(string(v)=>v for v in refs(i.expression,:param))
    isequal(expected,i.variables) || throw(ArgumentError("interaction $(i.name) variables do not match its expression"))
    all(k -> haskey(expected,k),keys(i.parameters)) ||
        throw(ArgumentError("interaction $(i.name) defaults contain absent variables"))
    declarations = Dict{Num,Number}(i.variables[k]=>v for (k,v) in i.parameters)
    validate_parameter_values(declarations)
    return declarations
end
