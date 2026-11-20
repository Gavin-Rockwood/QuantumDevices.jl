# Immediate upstream entries preserve declared tensor order without resolving names.
recipe_inputs(::Leaf) = NamedTuple[]
recipe_inputs(recipe::CompositeRecipe) = [
    (; name, recipe = binding.recipe, states = get(recipe.states, name, nothing))
    for (name, binding) in pairs(recipe.inputs)
]
recipe_inputs(recipe::Combine) = [
    (; name = input.name, recipe = input.recipe, states = input.states)
    for input in recipe.inputs
]
function recipe_inputs(recipe::ResolvedRecipe)
    selections = Dict(input.name => input.states for input in recipe_inputs(recipe.definition))
    return [(; name, recipe = child, states = get(selections, name, nothing))
            for (name, child) in recipe.inputs]
end

# Gather owner and operator declarations together, detecting only active-path cycles.
function collect_input_declarations(device, source, active = Any[])
    if source isa Symbol
        isnothing(device) && throw(ArgumentError("supply the device to inherit named recipes"))
        source = lookup(device.recipes, source, "recipe")
    end
    source isa Recipe || throw(ArgumentError("input must be a Recipe or registered name"))
    any(item -> item === source, active) && throw(ArgumentError("recipe cycle in input declarations"))
    if source isa Leaf
        names = if isnothing(device)
            nothing
        else
            component = lookup(device.components, source.component, "component")
            Symbol.(collect(keys(component.operators)))
        end
        return (; owners = [source.component], operators = Dict(source.component => names))
    end
    owners = Symbol[]
    operators = Dict{Symbol, Any}()
    for input in recipe_inputs(source)
        declarations = if source isa CompositeRecipe && isnothing(device) && input.recipe isa Symbol
            binding = source.inputs[input.name]
            (; owners = collect(Symbol, binding.owners), operators = Dict(binding.operator_names))
        else
            collect_input_declarations(device, input.recipe, vcat(active, Any[source]))
        end
        isempty(intersect(owners, declarations.owners)) ||
            throw(ArgumentError("overlapping input components at $(input.name)"))
        append!(owners, declarations.owners)
        merge!(operators, declarations.operators)
        if source isa CompositeRecipe && isnothing(device)
            binding = source.inputs[input.name]
            for (owner, names) in binding.operator_names
                operators[owner] = collect(Symbol, names)
            end
        end
    end
    return (; owners, operators)
end

# Owner-only legacy consumers can inspect leaves without component declarations.
recipe_owners(device, source) = collect_input_declarations(device, source).owners

# Component objects must be the device's actual registered declaration.

function component_name(device, component)
    name = component isa Component ? component.name : component
    name isa Symbol || throw(ArgumentError("component must be a registered Component or Symbol"))
    registered = lookup(device.components, name, "component")
    component isa Component &&
        component !== registered &&
        throw(ArgumentError("component object is not registered"))
    return name
end

# Normalize immutable-in-meaning selection requests, leaving occupation lookup to build time.
function retained_states(states)
    states isa AbstractVector || states isa AbstractRange ||
        throw(ArgumentError("states must be an ordered vector of zero-based indices or occupation tuples"))
    isempty(states) && throw(ArgumentError("states must be nonempty"))
    if all(i -> i isa Integer && !(i isa Bool), states)
        result = Int[states...]
        allunique(result) && all(>=(0), result) ||
            throw(ArgumentError("states must be distinct and nonnegative"))
        return result
    end
    all(i -> i isa Tuple, states) ||
        throw(ArgumentError("states cannot mix energy indices and occupation tuples"))
    width = length(first(states))
    width > 0 && all(t -> length(t) == width &&
        all(i -> i isa Integer && !(i isa Bool) && i >= 0, t), states) ||
        throw(ArgumentError("occupation tuples must have consistent nonzero length and nonnegative integer levels"))
    result = [Tuple(Int(i) for i in t) for t in states]
    allunique(result) || throw(ArgumentError("occupation tuples must be distinct"))
    return result
end


# Shared by execution and inspection. A single already-built recipe input with
# no new interaction is projected; automatic leaf initialization keeps its
# existing diagonalization/labeling path. Named inputs need device resolution.
function recipe_stage_operation(definition, sources)
    definition isa Leaf && return :initialize
    length(sources) == 1 && isempty(definition.interactions) || return :combine
    upstream = only(sources)
    upstream isa ResolvedRecipe && (upstream = upstream.definition)
    upstream isa Symbol && return :unresolved
    return upstream isa Leaf ? :combine : :projection
end
