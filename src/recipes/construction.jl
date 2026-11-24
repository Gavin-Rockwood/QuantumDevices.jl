"""
    Recipe(device, components; inputs=(;), interactions=Interaction[], initial_dims=Dict(), states=Dict(),
           parameters=Dict(), max_dimension=1000, dressing=true, dressing_options=(;))

Define a model from registered components without doing numerical work. Explicit
`inputs` aliases come first in tensor order, followed by uncovered components in
`components` order. Uncovered components are diagonalized at their default
dimension, or the dimension in `initial_dims`.

`states` selects ordered, zero-based eigenstates or dressed occupation tuples of each input, including grouped
inputs; omitted selections retain all states. Use component keys for automatic
inputs and alias symbols for explicit inputs. `interactions` accepts a collection
of `Interaction` objects using physical component operators, for example
`Recipe(device, [q1, q2]; interactions=[Interaction(:coupling,
0.02*q1.operators["x"]*q2.operators["x"])])`. Operators are bound to their
owning immediate input, including grouped inputs. Registered device interactions
can also be selected with `interactions=:auto`, for example
`Recipe(device, [q1, coupler, q2]; interactions=:auto)`. This selects registered
interactions involving at least two components, with every participant selected,
and spanning at least two immediate inputs. Terms internal to a grouped input
belong in that upstream recipe. Selection is captured at construction in sorted
interaction-name order; later device registrations do not alter it. The default
includes no interactions, and registered names are not accepted as a collection.
The supplied interactions are copied. To add interactions afterward,
use `add_interaction!(recipe, interaction)` with local interactions referencing
`recipe.inputs.alias.operators` (and then the component name for grouped inputs).
Only recipes registered with `add_recipe!` need repository names. Registration
captures an immutable upstream tree. Use `draft(saved_recipe)` to edit a copy,
and `replace_recipe!` to explicitly replace a registered definition.
`max_dimension=0` disables the tensor size limit. By default, physical occupation
labels follow a 20-sample cubic ramp of this node's new interactions. Set
`dressing=false` to diagonalize only the final Hamiltonian without tracking or
occupation labels. `dressing_options=(;)` configures enabled dressing
and is stored and validated even when dressing is disabled. Selections in
`states` supplied as integers remain zero-based energy indices. A vector of tuples
selects physical occupation labels in upstream component order, resolved afresh
after parameter overrides on each build. Missing labels require a larger upstream
retained space; label selection requires upstream dressing.

`Recipe(upstream; states=selection)` inherits upstream component and operator
declarations and projects its dressed eigenstates without repeating its internal
interactions. `Recipe((left=recipe1, right=recipe2); states=Dict(:left=>selection),
interactions=[...])` combines disjoint upstream recipes. Multiple inputs require
alias-keyed selections. Supply the device as the first argument for registered
recipe names or upstream declarations needing component lookup. For example,
`Recipe(intermediate; states=[(i,j,k) for i in 0:4 for j in 0:4 for k in 0:4],
max_dimension=125)` builds the full intermediate before selecting 125 labels.

`Recipe` is also the abstract supertype of all build definitions. `Leaf`, `Combine`,
and `Input` are lower-level compatibility interfaces.
"""
# Public construction normalizes inputs before shared interaction/default finalization.
Recipe(device, components; kwargs...) = CompositeRecipe(device, components; kwargs...)

# Validate stored dressing options even when dressing is disabled.
function validate_recipe_dressing(dressing, options)
    dressing isa Bool || throw(ArgumentError("dressing must be a Bool"))
    return normalize_dressing_options(options)
end

# Inherited construction captures declarations, never numerical state indices.
function inherited_recipe(device, inputs::NamedTuple;
                          states = Dict(), interactions = Interaction[],
                          parameters = Dict(), max_dimension::Int = 1000,
                          dressing = true, dressing_options = (;), control_maps=:auto)
    isempty(inputs) && throw(ArgumentError("inputs must be nonempty"))
    max_dimension >= 0 || throw(ArgumentError("max_dimension must be nonnegative"))
    scope = gensym(:recipe)
    owners = Symbol[]
    bindings = InputBinding[]
    for (alias, source) in pairs(inputs)
        source isa Union{Recipe, Symbol} || throw(ArgumentError("input $alias must be a recipe or registered name"))
        declarations = collect_input_declarations(device, source)
        operators = declarations.operators
        any(isnothing, values(operators)) && throw(ArgumentError(
            "input $alias: supply the device to inherit component operator declarations from $(typeof(source))"))
        localowners = declarations.owners
        allunique(localowners) || throw(ArgumentError("overlapping upstream input components at $alias"))
        isempty(intersect(owners, localowners)) || throw(ArgumentError("overlapping input components at $alias"))
        append!(owners, localowners)
        push!(bindings, InputBinding(alias, isnothing(device) ? source : copy_recipe_tree(device,source,false), localowners, scope, operators))
    end
    selections = Dict{Symbol, Any}()
    requests = if states isa AbstractDict
        states
    elseif length(inputs) == 1
        Dict(only(keys(inputs)) => states)
    else
        throw(ArgumentError("multiple inputs require states keyed by alias"))
    end
    for (alias, selection) in requests
        alias in keys(inputs) || throw(ArgumentError("unknown states input $alias"))
        selections[alias] = retained_states(selection)
    end
    return finish_recipe(device, owners, bindings, selections;
        parameters, interactions, max_dimension, scope, dressing, dressing_options, control_maps)
end

# Public inherited forms share finalization with explicit component construction.
Recipe(source::Recipe; kwargs...) = inherited_recipe(nothing, (; upstream = source); kwargs...)
Recipe(inputs::NamedTuple; kwargs...) = inherited_recipe(nothing, inputs; kwargs...)
Recipe(device, source::Union{Recipe, Symbol}; kwargs...) = inherited_recipe(device, (; upstream = source); kwargs...)
Recipe(device, inputs::NamedTuple; kwargs...) = inherited_recipe(device, inputs; kwargs...)

# Explicit aliases precede uncovered component leaves in tensor order.
function CompositeRecipe(
    device,
    components;
    inputs = (;),
    interactions = Interaction[],
    initial_dims = Dict(),
    states = Dict(),
    parameters = Dict(),
    max_dimension::Int = 1000,
    dressing = true,
    dressing_options = (;),
    control_maps = :auto,
)
    max_dimension >= 0 || throw(ArgumentError("max_dimension must be nonnegative"))
    inputs isa NamedTuple ||
        throw(ArgumentError("inputs must be a NamedTuple of aliases to recipes"))
    names = [component_name(device, c) for c in components]
    !isempty(names) && allunique(names) ||
        throw(ArgumentError("components must be nonempty and distinct"))
    scope = gensym(:recipe)
    bindings = InputBinding[]
    covered = Symbol[]
    for (alias, source) in pairs(inputs)
        source isa Union{Recipe, Symbol} ||
            throw(ArgumentError("input $alias must be a Recipe or registered name"))
        declarations = collect_input_declarations(device, source)
        owners = declarations.owners
        issubset(Set(owners), Set(names)) ||
            throw(ArgumentError("input $alias owns components outside this composite"))
        isempty(intersect(covered, owners)) ||
            throw(ArgumentError("overlapping input components at $alias"))
        append!(covered, owners)
        operators = declarations.operators
        push!(bindings, InputBinding(alias, copy_recipe_tree(device,source,false), owners, scope, operators))
    end
    automatic = [n for n in names if !(n in covered)]
    isempty(intersect(collect(keys(inputs)), automatic)) ||
        throw(ArgumentError("input alias collides with an automatic component input"))
    dimensions = Dict{Symbol, Int}()
    for (key, value) in initial_dims
        name = component_name(device, key)
        name in automatic ||
            throw(ArgumentError("initial_dims applies only to automatic leaves: $name"))
        value isa Integer && !(value isa Bool) && value > 0 ||
            throw(ArgumentError("initial dimension must be a positive integer"))
        haskey(dimensions, name) && throw(ArgumentError("duplicate initial dimension for $name"))
        dimensions[name] = value
    end
    for name in automatic
        component = lookup(device.components, name, "component")
        operators = Dict(name => Symbol.(collect(keys(component.operators))))
        push!(
            bindings,
            InputBinding(
                name,
                Leaf(component = name, dimension = get(dimensions, name, nothing), input_defaults=capture_leaf_inputs(component; defaults=device.defaults)),
                [name],
                scope,
                operators,
            ),
        )
    end
    selections = Dict{Symbol, Any}()
    aliases = [b.name for b in bindings]
    for (key, indices) in states
        alias = key isa Component ? component_name(device, key) : key
        alias isa Symbol && alias in aliases || throw(ArgumentError("unknown states input $alias"))
        key isa Component &&
            !(alias in automatic) &&
            throw(ArgumentError("component states apply only to automatic leaves"))
        haskey(selections, alias) && throw(ArgumentError("duplicate states input $alias"))
        selections[alias] = retained_states(indices)
    end
    return finish_recipe(device, names, bindings, selections;
        parameters, interactions, max_dimension, scope, dressing, dressing_options, control_maps)
end

# Both public paths capture the same ordered bindings, values, interactions, and available defaults.
function finish_recipe(device, components, bindings, states;
                       parameters, interactions, max_dimension, scope, dressing, dressing_options, control_maps)
    aliases = Tuple(binding.name for binding in bindings)
    recipe = CompositeRecipe(components, NamedTuple{aliases}(Tuple(bindings)), states,
        Dict{Num, Number}(parameters), Interaction[], max_dimension, scope, dressing, dressing_options,
        _capture_recipe_maps(device,bindings,control_maps,components))
    attach_interactions!(device, recipe, interactions)
    isnothing(device) || validate_component_input_defaults(device.components)
    return recipe
end
