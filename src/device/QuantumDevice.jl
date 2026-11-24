"""Named components, interactions, reusable build recipes, and nested gate variants."""
struct QuantumDevice{G<:AbstractDict}
    name::Symbol
    components::Dict{Symbol, Component}
    interactions::Dict{Symbol, Interaction}
    recipes::Dict{Symbol, Recipe}
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
        Dict{Symbol,Recipe}(), registry,
        InteractionMembership(:components,registry), InteractionMembership(:interactions,registry), metadata,
        Dict{Symbol,AbstractControlMap}(), Dict{Symbol,Union{Dict,QuantumGate}}())
end

# Reconstruct derived indexes through the same validated registration path.
function QuantumDevice(name, components::AbstractDict, interactions::AbstractDict,
    recipes::AbstractDict, metadata::NamedTuple)
    device = QuantumDevice(name; metadata)
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
    haskey(context,:explicit_parameters) && !isempty(context.explicit_parameters) &&
        throw(ArgumentError("old context defaults require recipe-value migration through load_device"))
    device = QuantumDevice(context.name; metadata=thaw_recipe_value(context.metadata))
    for c in values(context.components)
        component = Component(c.name; variables=c.variables, parameters=c.parameters,
            operators=c.operators, operator_generators=c.operator_generators,
            hamiltonian=c.hamiltonian, dimension=c.dimension, type=c.type, label=c.label,
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

# Removed default repositories must never be silently discarded by an old caller.
function QuantumDevice(name,components::AbstractDict,interactions::AbstractDict,recipes::AbstractDict,defaults::AbstractDict,metadata::NamedTuple)
    isempty(defaults) || throw(ArgumentError("device defaults were removed; supply values in owning recipe parameters"))
    return QuantumDevice(name,components,interactions,recipes,metadata)
end
