# Registered names and captured interaction declarations share display/default lookup.
interaction_name(interaction::Interaction) = interaction.name
interaction_name(name::Symbol) = name

# A grouped input owns its internal Hamiltonian. Only add couplings across
# immediate inputs here, using the device's captured full component membership.
function automatic_interactions(device, recipe::CompositeRecipe, selection::Symbol)
    selection === :auto || throw(ArgumentError(
        "unknown interactions setting: $selection; use :auto or a collection of Interaction objects"))
    input_for = Dict(owner => binding.name for binding in values(recipe.inputs) for owner in binding.owners)
    selected = Interaction[]
    for name in sort!(collect(keys(device.interaction_components)); by = string)
        owners = device.interaction_components[name]
        length(owners) >= 2 || continue
        all(owner -> haskey(input_for, owner), owners) || continue
        length(Set(input_for[owner] for owner in owners)) >= 2 || continue
        push!(selected, lookup(device.interactions, name, "interaction"))
    end
    return selected
end

# Preserve the symbolic tree, including ordered products and adjoints, while
# assigning each physical operator to its uniquely owning immediate input.
function bind_interaction_expression(expression::Sym, recipe::CompositeRecipe)
    if expression.expr[1] == :op
        ref = expression.expr[2]
        if ref isa OperatorRef
            owners = [binding for binding in values(recipe.inputs) if ref.owner in binding.owners]
            length(owners) == 1 ||
                throw(ArgumentError("interaction operator $(ref.owner).$(ref.name) must belong to exactly one recipe input"))
            binding = only(owners)
            ref.name in binding.operator_names[ref.owner] ||
                throw(ArgumentError("unknown input operator $(ref.owner).$(ref.name)"))
            return Sym((:op, BoundOperatorRef(recipe.scope, binding.name, ref)))
        end
        # add_interaction! validates existing bound refs without stripping scope.
        return expression
    end
    return Sym(Tuple(item isa Sym ? bind_interaction_expression(item, recipe) : item
        for item in expression.expr))
end

"""Add a coupling scoped exclusively to this composite's symbolic input operators."""
function add_interaction!(recipe::CompositeRecipe, interaction::Interaction)
    any(i -> i.name == interaction.name, recipe.interactions) &&
        throw(ArgumentError("duplicate local interaction name $(interaction.name)"))
    for ref in refs(interaction.expression, :op)
        ref isa BoundOperatorRef && ref.scope == recipe.scope && haskey(recipe.inputs, ref.input) ||
            throw(ArgumentError("local interactions must reference this recipe's input operators"))
        binding = recipe.inputs[ref.input]
        ref.physical.owner in binding.owners &&
        ref.physical.name in binding.operator_names[ref.physical.owner] ||
            throw(ArgumentError("unresolved input operator"))
    end
    validate_interaction_inputs(interaction)
    push!(recipe.interactions, interaction)
    return interaction
end


# Copy and bind each supplied interaction once into this recipe's local scope.
function attach_interactions!(device, recipe, interactions)
    selected = if interactions isa Symbol
        isnothing(device) && throw(ArgumentError("supply the device for automatic interactions"))
        automatic_interactions(device, recipe, interactions)
    else
        interactions
    end
    for interaction in selected
        interaction isa Interaction || throw(ArgumentError("interactions must contain Interaction objects"))
        captured = copy_interaction(interaction, false)
        expression = bind_interaction_expression(captured.expression, recipe)
        add_interaction!(recipe, Interaction(captured.name, expression, captured.variables, captured.defaults, captured.metadata))
    end
    return recipe
end

