"""Register a component after validating owners, references, and parameters."""
function add_component!(device, component)
    name = component.name
    (haskey(device.components, name) || haskey(device.interactions, name)) &&
        throw(ArgumentError("name already used"))
    # Revalidate mutable public declaration maps at registration.
    Component(name; variables=component.variables, parameters=component.parameters,
        operators=component.operators, operator_generators=component.operator_generators,
        hamiltonian=component.hamiltonian, dimension=component.dimension, type=component.type, label=component.label,
        metadata=component.metadata)
    validate_device_defaults(merge(device.components,Dict(component.name=>component)),device.interactions)
    physical = [physical_operator(operator) for operator in values(component.operators)]
    for ref in physical
        haskey(device.operators.by_operator,ref) &&
            throw(ArgumentError("operator $ref already belongs to $(device.operators.by_operator[ref])"))
    end
    component.operators.entries isa Tuple || freeze_operator_declarations!(component.operators)
    component.operator_generators.entries isa Tuple || freeze_operator_declarations!(component.operator_generators)
    device.components[name] = component
    getfield(device.operators,:registered_components)[name] = component
    getfield(device.operators,:by_component)[name] = component.operators
    for ref in physical
        getfield(device.operators,:by_operator)[ref] = name
    end
    getfield(device.interactions_for,:entries)[name] = ()
    return component
end

"""Register an interaction whose operators belong to registered components."""
function add_interaction!(device, interaction)
    name = interaction.name
    (haskey(device.components, name) || haskey(device.interactions, name)) &&
        throw(ArgumentError("name already used"))
    validate_interaction_inputs(interaction)
    all(ref -> ref isa OperatorRef, refs(interaction.expression,:op)) ||
        throw(ArgumentError("device interactions require physical component operator references"))
    validate_device_defaults(device.components,merge(device.interactions,Dict(name=>interaction)))
    owners = Tuple(get_components(interaction,device))
    # Membership depends only on the immutable symbolic expression. Capture declaration
    # variables and numerical declaration defaults as well.
    captured = Interaction(interaction.name, interaction.expression,
        RecipeMap(interaction.variables), Dict{Num,Union{Missing,Number}}(interaction.defaults), interaction.metadata)
    device.interactions[name] = captured
    getfield(device.interaction_components,:entries)[name] = owners
    for owner in owners
        getfield(device.interactions_for,:entries)[owner] = (device.interactions_for[owner]...,name)
    end
    return captured
end

"""Register a uniquely named reusable recipe."""
function add_recipe!(device, name, recipe::Recipe)
    haskey(device.recipes, name) && throw(ArgumentError("recipe already exists"))
    validate_component_input_defaults(device.components)
    components = Dict{Symbol,Component}()
    captured = copy_recipe_tree(device,recipe,true;components)
    resolve_recipe_inputs(captured,components)
    device.recipes[name] = captured
    return captured
end

"""Replace an existing registered recipe with a freshly captured immutable tree."""
function replace_recipe!(device, name, recipe::Recipe)
    haskey(device.recipes, name) || throw(ArgumentError("unknown recipe: $name"))
    validate_component_input_defaults(device.components)
    components = Dict{Symbol,Component}()
    captured = copy_recipe_tree(device,recipe,true;components)
    resolve_recipe_inputs(captured,components)
    device.recipes[name] = captured
    return captured
end

# A parameter symbol identifies one physical value throughout a device.
function validate_parameter_defaults(parameters)
    seen = Dict{Num,Number}()
    for (variable, value) in parameters
        is_parameter(variable) && concrete_parameter_value(value) ||
            throw(ArgumentError("parameter declarations require atomic variables and finite concrete numbers"))
        haskey(seen, variable) && !isequal(seen[variable], value) &&
            throw(ArgumentError("conflicting defaults for parameter $variable: $(seen[variable]) and $value"))
        seen[variable] = value
    end
    return nothing
end

function lookup(repository::AbstractDict{Symbol}, name::Symbol, kind)
    haskey(repository, name) || throw(ArgumentError("unknown $kind: $name"))
    return repository[name]
end

# Validate declarations before applying any operating-point overrides.
function validate_default_sources(groups)
    seen = Dict{Num,Tuple{Number,String}}()
    for (source, declarations) in groups, (variable,value) in declarations
        validate_parameter_values(Dict(variable=>value))
        if haskey(seen,variable)
            previous, origin = seen[variable]
            isequal(previous,value) || throw(ArgumentError(
                "conflicting defaults for parameter $variable: $previous from $origin and $value from $source"))
        else
            seen[variable] = (value,string(source))
        end
    end
    return nothing
end

# Numerical defaults are derived from declarations, independently of recipe usage.
function validate_device_defaults(components,interactions)
    groups = Any[("component $(c.name)",component_parameter_values(c)) for c in values(components)]
    for term in values(interactions)
        validate_interaction_inputs(term)
        push!(groups,("interaction $(term.name)",Dict(v=>x for (v,x) in term.defaults if !ismissing(x))))
    end
    validate_default_sources(groups)
    return nothing
end
function device_defaults(device)
    validate_device_defaults(device.components,device.interactions)
    result = Dict{String,Union{Missing,Number}}()
    for component in values(device.components)
        idle = component_parameter_values(component)
        for variable in values(component.variables)
            name = string(variable)
            value = get(idle,variable,missing)
            (!haskey(result,name) || ismissing(result[name])) && (result[name]=value)
        end
    end
    for term in values(device.interactions), (variable,value) in term.defaults
        name = string(variable)
        (!haskey(result,name) || ismissing(result[name])) && (result[name]=value)
    end
    return result
end
