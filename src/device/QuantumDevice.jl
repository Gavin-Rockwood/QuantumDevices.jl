"""Named components, interactions, reusable build recipes, and nested gate variants."""
struct QuantumDevice{G<:AbstractDict}
    name::Symbol
    components::Dict{Symbol, Component}
    interactions::Dict{Symbol, Interaction}
    recipes::Dict{Symbol, Recipe}
    explicit_parameters::Dict{Num,Number}
    operators::OperatorRegistry
    interactions_for::InteractionMembership
    interaction_components::InteractionMembership
    metadata::NamedTuple
    control_maps::Dict{Symbol,AbstractControlMap}
    gates::G
end

function QuantumDevice(name; metadata = (;))
    registry = OperatorRegistry()
    return QuantumDevice(name, Dict{Symbol,Component}(), Dict{Symbol,Interaction}(),
        Dict{Symbol,Recipe}(), Dict{Num,Number}(), registry,
        InteractionMembership(:components,registry), InteractionMembership(:interactions,registry), metadata,
        Dict{Symbol,AbstractControlMap}(), Dict{Symbol,Union{Dict,QuantumGate}}())
end

# Reconstruct derived indexes through the same validated registration path.
function QuantumDevice(name, components::AbstractDict, interactions::AbstractDict,
    recipes::AbstractDict, explicit_parameters::AbstractDict, metadata::NamedTuple)
    device = QuantumDevice(name; metadata)
    add_parameters!(device,explicit_parameters)
    for component in values(components)
        add_component!(device,component)
    end
    for interaction in values(interactions)
        add_interaction!(device,interaction)
    end
    merge!(device.recipes,recipes)
    return device
end

"""Rebuild an independent device and its indexes from a saved physical gate context."""
function QuantumDevice(context::NamedTuple)
    device = QuantumDevice(context.name; metadata=thaw_recipe_value(context.metadata))
    add_parameters!(device,Dict(context.explicit_parameters))
    for c in values(context.components)
        component = Component(c.name; variables=c.variables, parameters=c.parameters,
            operators=c.operators, operator_generators=c.operator_generators,
            hamiltonian=c.hamiltonian, type=c.type, label=c.label,
            metadata=thaw_recipe_value(c.metadata))
        add_component!(device,component)
    end
    for interaction in values(context.interactions)
        add_interaction!(device,copy_interaction(interaction,false))
    end
    for (name,map) in get(context,:control_maps,Dict())
        name === map.name || throw(ArgumentError("control map registry name does not match $name"))
        add_control_map!(device,map)
    end
    merge!(device.recipes,context.recipes)
    return device
end

# Preserve the previous direct field constructor; derived indexes remain supplied by its caller.
QuantumDevice(name::Symbol,components::Dict{Symbol,Component},interactions::Dict{Symbol,Interaction},
    recipes::Dict{Symbol,Recipe},explicit_parameters::Dict{Num,Number},operators::OperatorRegistry,
    interactions_for::InteractionMembership,interaction_components::InteractionMembership,metadata::NamedTuple) =
    QuantumDevice(name,components,interactions,recipes,explicit_parameters,operators,
        interactions_for,interaction_components,metadata,Dict{Symbol,AbstractControlMap}(),Dict{Symbol,Union{Dict,QuantumGate}}())

# Preserve the previous complete field constructor with an empty registry.
QuantumDevice(name,components,interactions,recipes,explicit_parameters,operators,interactions_for,interaction_components,metadata,gates::AbstractDict) =
    QuantumDevice(name,components,interactions,recipes,explicit_parameters,operators,interactions_for,interaction_components,metadata,Dict{Symbol,AbstractControlMap}(),gates)
