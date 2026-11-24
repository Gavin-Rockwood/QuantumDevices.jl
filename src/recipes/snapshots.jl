# Rebind only operators captured in the copied node's scope.
function copy_expression(expression::Sym, old_scope, new_scope)
    if expression.expr[1] == :op
        ref = expression.expr[2]
        if ref isa BoundOperatorRef && ref.scope == old_scope
            return Sym((:op, BoundOperatorRef(new_scope, ref.input, ref.physical)))
        end
    end
    return Sym(Tuple(item isa Sym ? copy_expression(item, old_scope, new_scope) : item
                     for item in expression.expr))
end

# Snapshot collections are frozen; editable drafts thaw the same captured values.
function copy_interaction(interaction::Interaction, frozen; old_scope = nothing, new_scope = nothing)
    metadata = freeze_recipe_value(interaction.metadata)
    if !frozen
        metadata = thaw_recipe_value(metadata)
    end
    variables = frozen ? RecipeMap(interaction.variables) : Dict{String, Num}(interaction.variables)
    expression = copy_expression(interaction.expression, old_scope, new_scope)
    return Interaction(interaction.name, expression, variables, metadata)
end

# Frozen storage is a provenance property, not proof of device compatibility.
is_recipe_snapshot(recipe::Recipe) = recipe.overrides isa RecipeMap

# Capture errors identify the occurrence being copied without a diagnostic object.
capture_error(path, message) = throw(ArgumentError(join(path, " -> ") * ": " * message))
function capture_lookup(repository, name, kind, path)
    haskey(repository, name) || capture_error(path, "unknown $kind: $name")
    return repository[name]
end

# Capture each current component declaration once per build, without generating matrices.
function capture_component!(components, device, name, path)
    return get!(components, name) do
        component = capture_lookup(device.components, name, "component", path)
        component.name == name || capture_error(path, "component registry name does not match $name")
        try
            Component(name;
                variables = component.variables,
                parameters = component.parameters,
                operators = component.operators,
                operator_generators = component.operator_generators,
                hamiltonian = component.hamiltonian,
                dimension = component.dimension,
                type = component.type,
                label = component.label,
                metadata = component.metadata)
        catch error
            error isa ArgumentError || rethrow()
            capture_error(path, error.msg)
        end
    end
end

# Local interactions require scoped operators; gate additions may use physical refs.
function capture_operator_refs(interaction, definition, owners, operators, path; local_scope)
    for ref in refs(interaction.expression, :op)
        ref isa Union{OperatorRef, BoundOperatorRef} || capture_error(path, "invalid operator reference")
        if ref isa BoundOperatorRef
            definition isa CompositeRecipe && ref.scope == definition.scope &&
                haskey(definition.inputs, ref.input) || capture_error(path, "foreign input operator scope")
            binding = definition.inputs[ref.input]
            ref.physical.owner in binding.owners &&
                ref.physical.name in binding.operator_names[ref.physical.owner] ||
                capture_error(path, "unresolved input operator")
        elseif local_scope && definition isa CompositeRecipe
            capture_error(path, "local interactions require scoped input operators")
        end
        physical = physical_ref(ref)
        physical.owner in owners || capture_error(path, "interaction outside selected inputs")
        current = operators[physical.owner]
        isnothing(current) || physical.name in current ||
            capture_error(path, "unknown input operator $(physical.owner).$(physical.name)")
    end
    return nothing
end

# Public/internal copy callers need only the recipe; capture consumers also use its declarations.
function copy_recipe_tree(device, source, frozen; components = Dict{Symbol, Component}())
    return capture_recipe_tree(device, source, frozen, Any[], ["target"], components).recipe
end

# Copy and structural validation share one contextual traversal. Child declarations
# and interaction provenance are temporary summaries, never persistent validity flags.
function capture_recipe_tree(device, source, frozen, active, path, components)
    definition = source
    if source isa Symbol
        isnothing(device) && capture_error(path, "supply the device to capture named recipe inputs")
        definition = capture_lookup(device.recipes, source, "recipe", path)
    end
    definition isa Recipe || capture_error(path, "expected a Recipe or registered name")
    any(item -> item === definition, active) && capture_error(path, "recipe cycle during capture")
    active = vcat(active, Any[definition])
    parameters = frozen ? RecipeMap(definition.overrides) : Dict{Num, Number}(definition.overrides)
    if definition isa Leaf
        operators = Dict{Symbol, Any}(definition.component => nothing)
        if !isnothing(device)
            component = capture_component!(components, device, definition.component, path)
            isnothing(definition.input_defaults) || Set(keys(definition.input_defaults)) == Set(values(component.variables)) ||
                capture_error(path,"component input declarations changed since recipe capture")
            operators[definition.component] = Set(Symbol.(keys(component.operators)))
        end
        for (ref, generator) in definition.generators
            ref isa OperatorRef && ref.owner == definition.component ||
                capture_error(path, "generator override must belong to $(definition.component)")
            generator isa Function || capture_error(path, "invalid operator generator override")
            applicable(generator, something(definition.dimension, 1)) ||
                capture_error(path, "operator generator must accept a dimension")
            names = operators[definition.component]
            isnothing(names) || ref.name in names ||
                capture_error(path, "override refers to unknown operator $(ref.name)")
        end
        generators = Dict{OperatorRef, Function}(
            ref => generator
            for (ref, generator) in definition.generators)
        captured = Leaf(definition.component, definition.dimension,
                        frozen ? RecipeMap(generators) : generators, parameters,
                        _copy_recipe_maps(device,definition.control_maps,[definition.component],frozen),
                        isnothing(definition.input_defaults) ? (isnothing(device) ? nothing : capture_leaf_inputs(components[definition.component])) : definition.input_defaults)
        return (; recipe = captured, owners = [definition.component], operators, provenance = Set{Any}())
    end
    definition isa Union{CompositeRecipe, Combine} ||
        capture_error(path, "unsupported recipe $(typeof(definition))")
    composite = definition isa CompositeRecipe
    bindings = composite ? collect(values(definition.inputs)) : collect(definition.inputs)
    all(input -> composite ? input isa InputBinding : input isa Input, bindings) ||
        capture_error(path, "invalid input binding declaration")
    entries = recipe_inputs(definition)
    aliases = [input.name for input in entries]
    isempty(entries) && capture_error(path, "inputs must be nonempty")
    allunique(aliases) || capture_error(path, "duplicate local input names")
    if composite
        allunique(definition.components) && !isempty(definition.components) ||
            capture_error(path, "components must be nonempty and distinct")
    end
    owners = Symbol[]
    operators = Dict{Symbol, Any}()
    provenance = Set{Any}()
    copied = Any[]
    scope = composite ? gensym(:recipe) : nothing
    for (entry, input) in zip(entries, bindings)
        alias = entry.name
        input.name == alias || capture_error(path, "input alias $alias differs from binding $(input.name)")
        childpath = vcat(path, [string(alias)])
        child = capture_recipe_tree(device, entry.recipe, frozen, active, childpath, components)
        isempty(intersect(owners, child.owners)) || capture_error(childpath, "overlapping input components")
        isempty(intersect(provenance, child.provenance)) || capture_error(childpath, "duplicate interaction provenance")
        union!(provenance, child.provenance)
        append!(owners, child.owners)
        if composite
            input.scope == definition.scope || capture_error(childpath, "foreign input binding scope")
            allunique(input.owners) && Set(input.owners) == Set(child.owners) ||
                capture_error(childpath, "input ownership changed")
            Set(keys(input.operator_names)) == Set(child.owners) ||
                capture_error(childpath, "input operator owners changed")
            for owner in child.owners
                names = input.operator_names[owner]
                all(name -> name isa Symbol, names) && allunique(names) ||
                    capture_error(childpath, "invalid input operator names")
                current = child.operators[owner]
                isnothing(current) || Set(names) == current ||
                    capture_error(childpath, "input operators changed for $owner")
                operators[owner] = Set(names)
            end
            own = frozen ? RecipeVector(input.owners) : collect(Symbol, input.owners)
            names = frozen ? RecipeMap(input.operator_names) :
                Dict(owner => collect(Symbol, names) for (owner, names) in input.operator_names)
            push!(copied, InputBinding(alias, child.recipe, own, scope, names))
        else
            merge!(operators, child.operators)
            states = retained_states(entry.states)
            push!(copied, Input(alias, child.recipe, frozen ? RecipeVector(states) : states))
        end
    end
    if composite
        Set(owners) == Set(definition.components) ||
            capture_error(path, "input ownership changed since construction")
    end
    states = Dict{Symbol, Any}()
    if composite
        for (alias, selection) in definition.states
            alias in aliases || capture_error(path, "unknown states input $alias")
            states[alias] = retained_states(selection)
        end
    end
    interactions = Interaction[]
    names = Symbol[]
    for (index, item) in enumerate(definition.interactions)
        interaction = item
        if item isa Symbol
            isnothing(device) && capture_error(path, "supply the device for named interactions")
            interaction = capture_lookup(device.interactions, item, "interaction", path)
        end
        interaction isa Interaction || capture_error(path, "interactions must contain Interaction objects")
        interaction.name in names && capture_error(path, "duplicate local interaction name $(interaction.name)")
        push!(names, interaction.name)
        origin = composite ? (definition.scope, index) : (:device_interaction, interaction.name)
        origin in provenance && capture_error(path, "duplicate interaction provenance")
        push!(provenance, origin)
        capture_operator_refs(interaction, definition, owners, operators, path; local_scope = true)
        push!(interactions, copy_interaction(interaction, frozen;
            old_scope = composite ? definition.scope : nothing, new_scope = scope))
    end
    if composite
        captured = CompositeRecipe(
            frozen ? RecipeVector(definition.components) : collect(Symbol, definition.components),
            NamedTuple{Tuple(aliases)}(Tuple(copied)), frozen ? RecipeMap(states) : states,
            parameters, frozen ? RecipeVector(interactions) : interactions,
            definition.max_dimension, scope, definition.dressing, definition.dressing_options,
            _copy_recipe_maps(device,definition.control_maps,owners,frozen))
    else
        captured = Combine(frozen ? RecipeVector(copied) : Input[copied...],
            frozen ? RecipeVector(interactions) : Union{Symbol, Interaction}[interactions...],
            parameters, definition.max_dimension, definition.dressing, definition.dressing_options,
            _copy_recipe_maps(device,definition.control_maps,owners,frozen))
    end
    return (; recipe = captured, owners, operators, provenance)
end

# Device-less provenance may reuse frozen data; device builds always recheck compatibility.
function snapshot_recipe(device, recipe; components = Dict{Symbol, Component}())
    isnothing(device) && recipe isa Recipe && is_recipe_snapshot(recipe) && return recipe
    return copy_recipe_tree(device, recipe, true; components)
end

"""
    draft(saved_recipe)

Return an independent editable recipe tree. Scoped input operators are rebound
consistently. Callables retain their identity and may capture mutable external
state, which is outside the snapshot boundary.
"""
draft(recipe::Recipe) = copy_recipe_tree(nothing, recipe, false)
