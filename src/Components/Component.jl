# Stable identities are declared beside builtin generator functions.
const BUILTIN_OPERATOR_FUNCTIONS = Dict{String,Function}()

"""Normalize a positive construction Hilbert-space state count."""
function validate_operator_dimension(dimension)
    dimension isa Integer && dimension > 0 && dimension <= typemax(Int) ||
        throw(ArgumentError("dimension must be a positive representable integer state count"))
    return Int(dimension)
end

"""A symbolic component with matching string-keyed declaration maps."""
struct Component
    name::Symbol
    variables::Dict{String,Num}
    parameters::Dict{String,Number}
    operators::OperatorDeclarations{Sym}
    operator_generators::OperatorDeclarations{Function}
    hamiltonian::Sym
    dimension::Union{Nothing,Int}
    type::String
    label::Union{Nothing,String}
    metadata::NamedTuple
    function Component(name; variables, parameters, operators, operator_generators,
        hamiltonian, dimension=nothing, type="custom", label=string(name), metadata=(;))
        name isa Symbol || throw(ArgumentError("component name must be a Symbol"))
        variables = named_component_map(variables, Num, "variables")
        parameters = named_component_map(parameters, Number, "parameters")
        operators = named_component_map(operators, Sym, "operators")
        operator_generators = named_component_map(operator_generators, Function, "operator_generators")
        Set(keys(variables)) == Set(keys(parameters)) ||
            throw(ArgumentError("variables and parameters must have matching local keys"))
        Set(keys(operators)) == Set(keys(operator_generators)) ||
            throw(ArgumentError("operators and operator_generators must have matching local keys"))
        dimension = isnothing(dimension) ? nothing : validate_operator_dimension(dimension)
        all(generator -> applicable(generator, something(dimension, 1)), values(operator_generators)) ||
            throw(ArgumentError("operator generators must be callable with a dimension"))
        all(is_parameter, values(variables)) || throw(ArgumentError("variables must be atomic Symbolics scalar variables"))
        allunique(values(variables)) || throw(ArgumentError("duplicate component parameter symbols"))
        all(concrete_parameter_value, values(parameters)) || throw(ArgumentError("parameter values must be finite concrete numbers"))
        for (key, expression) in operators
            expression.expr == (:op, OperatorRef(name, Symbol(key))) ||
                throw(ArgumentError("operator $key must be the atomic operator owned by $name with that local name"))
        end
        hamiltonian isa Sym || throw(ArgumentError("hamiltonian must be a symbolic operator expression"))
        declared = Set(expression.expr[2] for expression in values(operators))
        issubset(refs(hamiltonian, :op), declared) || throw(ArgumentError("Hamiltonian contains undeclared or foreign operators"))
        issubset(refs(hamiltonian, :param), Set(values(variables))) || throw(ArgumentError("Hamiltonian parameters must be declared in variables"))
        return new(name, variables, parameters, OperatorDeclarations(operators), OperatorDeclarations(operator_generators),
            hamiltonian, dimension, string(type), isnothing(label) ? nothing : string(label), metadata)
    end
end

"""Normalize local Symbol/string keys, rejecting aliases that collide after normalization."""
function named_component_map(declarations, ::Type{T}, description) where T
    result = Dict{String,T}()
    for (key, value) in pairs(declarations)
        key isa Union{Symbol,AbstractString} || throw(ArgumentError("$description keys must be strings or Symbols"))
        local_name = string(key)
        isempty(local_name) && throw(ArgumentError("$description keys must be nonempty"))
        haskey(result, local_name) && throw(ArgumentError("duplicate $description key: $local_name"))
        value isa T || throw(ArgumentError("invalid $description value for $local_name"))
        result[local_name] = value
    end
    return result
end

"""Return symbolic operators and their factories as matching string-keyed maps."""
function component_operators(name::Symbol, declarations=NamedTuple(); kwargs...)
    generators = named_component_map(declarations, Function, "operator_generators")
    for (key, generator) in named_component_map(kwargs, Function, "operator_generators")
        haskey(generators, key) && throw(ArgumentError("duplicate operator_generators key: $key"))
        generators[key] = generator
    end
    return Dict(key => op(name, Symbol(key)) for key in keys(generators)), generators
end

# Resolve public local names only at the boundary where physical symbols are needed.
function component_parameter_values(c)
    Set(keys(c.variables)) == Set(keys(c.parameters)) ||
        throw(ArgumentError("component $(c.name) variables and parameters must have matching local keys"))
    allunique(values(c.variables)) || throw(ArgumentError("duplicate component parameter symbols"))
    result = Dict{Num,Number}(c.variables[key] => value for (key,value) in c.parameters)
    validate_parameter_values(result)
    return result
end
component_operator_generators(c) = Dict{OperatorRef,Function}(c.operators[key].expr[2] => generator for (key,generator) in c.operator_generators)

# Positional declarations remain accepted for custom constructors written before
# the named operator maps were introduced.
function Component(name, variables, parameters, generators, hamiltonian, type, label, metadata; dimension=nothing)
    vars = named_component_map(variables, Num, "variables")
    defaults = all(key -> key isa Num, keys(parameters)) ?
        Dict(key => parameters[value] for (key,value) in vars) : parameters
    ops = Dict(string(ref.name) => op(ref.owner, ref.name) for ref in keys(generators))
    gens = Dict(string(ref.name) => generator for (ref,generator) in generators)
    return Component(name; variables=vars, parameters=defaults, operators=ops,
        operator_generators=gens, hamiltonian, dimension, type, label, metadata)
end
