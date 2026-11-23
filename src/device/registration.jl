"""Register a component after validating owners, references, and parameters."""
function add_component!(device, component)
    name = component.name
    (haskey(device.components, name) || haskey(device.interactions, name)) &&
        throw(ArgumentError("name already used"))
    # Revalidate mutable public declaration maps at registration.
    Component(name; variables=component.variables, parameters=component.parameters,
        operators=component.operators, operator_generators=component.operator_generators,
        hamiltonian=component.hamiltonian, type=component.type, label=component.label,
        metadata=component.metadata)
    validate_default_sources(vcat(device_default_sources(device), [("component $(component.name)",component_parameter_values(component))]))
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
    validate_default_sources(vcat(device_default_sources(device), [("registered interaction $(interaction.name)",interaction_parameter_values(interaction))]))
    all(ref -> ref isa OperatorRef, refs(interaction.expression,:op)) ||
        throw(ArgumentError("device interactions require physical component operator references"))
    owners = Tuple(get_components(interaction,device))
    # Membership depends only on the immutable symbolic expression. Capture declaration
    # variables as well, while leaving numerical defaults editable as before.
    captured = Interaction(interaction.name, interaction.expression,
        RecipeMap(interaction.variables), interaction.parameters, interaction.metadata)
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
    components = Dict{Symbol,Component}()
    captured = copy_recipe_tree(device,recipe,true;components)
    node = prepare_captured_node(captured,components;enforce_limits=false)
    validate_default_sources(vcat(device_default_sources(device),node_default_sources(node)))
    device.recipes[name] = captured
    return captured
end

"""Replace an existing registered recipe with a freshly captured immutable tree."""
function replace_recipe!(device, name, recipe::Recipe)
    haskey(device.recipes, name) || throw(ArgumentError("unknown recipe: $name"))
    components = Dict{Symbol,Component}()
    captured = copy_recipe_tree(device,recipe,true;components)
    node = prepare_captured_node(captured,components;enforce_limits=false)
    validate_default_sources(vcat(device_default_sources(device),node_default_sources(node)))
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

function device_parameters(device)
    parameters = collect(pairs(device.explicit_parameters))
    for c in values(device.components)
        append!(parameters, pairs(component_parameter_values(c)))
    end
    for i in values(device.interactions)
        append!(parameters, pairs(interaction_parameter_values(i)))
    end
    return parameters
end

"""Assemble full symbolic names without discarding unassigned variables."""
function device_variables(device)
    variables = Dict{String,Num}()
    declarations = Iterators.flatten((
        (values(c.variables) for c in values(device.components)),
        (values(i.variables) for i in values(device.interactions)),
        (keys(device.explicit_parameters),),
    ))
    for group in declarations, variable in group
        name = string(variable)
        haskey(variables, name) && !isequal(variables[name], variable) &&
            throw(ArgumentError("distinct device variables have the same name: $name"))
        variables[name] = variable
    end
    return variables
end

# Pair every named device variable with its concrete default or `missing`.
function device_parameter_view(device)
    validate_default_sources(device_default_sources(device))
    defaults = Dict{Num,Number}(device_parameters(device))
    return Dict{String,Union{Missing,Number}}(
        name => get(defaults, variable, missing) for (name, variable) in device_variables(device))
end


"""Atomically add explicit shared defaults. Equal existing defaults are allowed.
`device.parameters` is an assembled view; use this function to declare defaults.
"""
function add_parameters!(device::QuantumDevice, parameters::AbstractDict)
    validate_parameter_values(parameters)
    candidate = merge(device.explicit_parameters, parameters)
    validate_default_sources(vcat(device_default_sources(device), [("new device explicit parameters",parameters)]))
    merge!(device.explicit_parameters, candidate)
    return device
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
function device_default_sources(device)
    groups = Any[("device $(device.name) explicit parameters", device.explicit_parameters)]
    append!(groups, [("component $(c.name)",component_parameter_values(c)) for c in values(device.components)])
    append!(groups, [("registered interaction $(i.name)",interaction_parameter_values(i)) for i in values(device.interactions)])
    return groups
end
