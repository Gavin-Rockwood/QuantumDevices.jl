"""
    QuantumDevice(name; components=(), interactions=(), modelspecs=(), gatespecs=())

Mutable symbolic registry for a physical device. Initial entries are registered
in dependency order: components, interactions, model recipes, then gate recipes.
"""
struct QuantumDevice
    name::Symbol
    components::Dict{Symbol,Component}
    interactions::Dict{Symbol,OperatorExpr}
    modelspecs::Dict{Symbol,ModelSpec}
    gatespecs::Dict{Symbol,GateSpec}
end

function QuantumDevice(
        name::Symbol;
        components = (),
        interactions = (),
        modelspecs = (),
        gatespecs = ())
    device = QuantumDevice(
        name,
        Dict{Symbol,Component}(),
        Dict{Symbol,OperatorExpr}(),
        Dict{Symbol,ModelSpec}(),
        Dict{Symbol,GateSpec}(),
    )

    for component in _device_contents(components, "components")
        component isa Component ||
            error("QuantumDevice components must be Component objects.")
        register!(device, component)
    end
    for entry in _device_contents(interactions, "interactions")
        entry isa Pair ||
            error("QuantumDevice interactions must be name => expression pairs.")
        first(entry) isa Symbol ||
            error("QuantumDevice interaction names must be symbols.")
        last(entry) isa OperatorExpr ||
            error("QuantumDevice interactions must be OperatorExpr objects.")
        register!(device, first(entry), last(entry))
    end
    for spec in _device_contents(modelspecs, "modelspecs")
        spec isa ModelSpec ||
            error("QuantumDevice modelspecs must be ModelSpec objects.")
        register!(device, spec)
    end
    for gate in _device_contents(gatespecs, "gatespecs")
        gate isa GateSpec ||
            error("QuantumDevice gatespecs must be GateSpec objects.")
        register!(device, gate)
    end
    device
end

function _device_contents(contents, registry)
    contents isa Union{Tuple,AbstractVector} ||
        error("QuantumDevice $registry must be a tuple or vector.")
    contents
end

function register!(device::QuantumDevice, component::Component)
    _require_new_name(device.components, component.name, "component")
    device.components[component.name] = component
    device
end

function register!(
        device::QuantumDevice,
        name::Symbol,
        interaction::OperatorExpr)
    _require_new_name(device.interactions, name, "interaction")
    device.interactions[name] = interaction
    device
end

function register!(device::QuantumDevice, spec::ModelSpec)
    _require_new_name(device.modelspecs, spec.name, "model")
    for component in values(spec.subsystems)
        _require_registered_component(device, component, spec)
    end
    for interaction in spec.interactions
        any(candidate -> candidate === interaction, values(device.interactions)) ||
            error("ModelSpec $(spec.name) contains an unregistered interaction expression.")
    end
    device.modelspecs[spec.name] = spec
    device
end

function register!(device::QuantumDevice, gate::GateSpec)
    _require_new_name(device.gatespecs, gate.name, "gate")
    haskey(device.modelspecs, gate.modelspec.name) &&
        device.modelspecs[gate.modelspec.name] === gate.modelspec ||
        error(
            "GateSpec $(gate.name) requires its exact ModelSpec " *
            "$(gate.modelspec.name) to be registered.",
        )
    device.gatespecs[gate.name] = gate
    device
end

function _require_new_name(registry, name, kind)
    haskey(registry, name) &&
        error("QuantumDevice already contains a $kind named $name.")
end

function _require_registered_component(device, component, spec)
    haskey(device.components, component.name) ||
        error("ModelSpec $(spec.name) contains unregistered component $(component.name).")
    registered = device.components[component.name]
    component === registered && return
    component.spec isa GenericSpec &&
        component.spec.source === registered.spec && return
    error(
        "ModelSpec $(spec.name) does not contain the registered " *
        "component object $(component.name).",
    )
end

"""
    modelspec(device, name, components, interactions=(); parameters=(), kwargs...)

Create an unsaved `ModelSpec` by resolving ordered component and interaction
names from `device`.
"""
function modelspec(
        device::QuantumDevice,
        name::Symbol,
        component_names,
        interaction_names = ();
        parameters = (),
        dims = Dict(),
        kwargs...)
    components = _resolve_device_names(
        device.components,
        component_names,
        "component",
    )
    interactions = _resolve_device_names(
        device.interactions,
        interaction_names,
        "interaction",
    )
    parameters isa Union{Tuple,AbstractVector} ||
        error("Model parameters must be a tuple or vector.")
    model_parameters = Dict{ParamPath,DeviceParameter}()
    for interaction in interactions
        for (path, parameter) in QuantumDevices.parameters(interaction)
            model_parameters[path] = parameter
        end
    end
    for parameter in parameters
        parameter isa DeviceParameter ||
            error("Model parameters must be DeviceParameter objects.")
        model_parameters[parameter.path] = parameter
    end
    ModelSpec(
        name,
        Component[components...],
        dims;
        interactions = Tuple(interactions),
        parameters = Tuple(values(model_parameters)),
        kwargs...,
    )
end

function _resolve_device_names(registry, names, kind)
    names isa Union{Tuple,AbstractVector} ||
        error("Model $kind names must be a tuple or vector.")
    resolved = Any[]
    for name in names
        name isa Symbol || error("Model $kind names must be symbols.")
        haskey(registry, name) ||
            error("QuantumDevice does not contain a $kind named $name.")
        push!(resolved, registry[name])
    end
    resolved
end

function model(device::QuantumDevice, name::Symbol)
    haskey(device.modelspecs, name) ||
        error("QuantumDevice does not contain a model named $name.")
    model(device.modelspecs[name])
end

function numerical(
        device::QuantumDevice,
        model_name::Symbol,
        expression::Union{OperatorExpr,ScalarExpr})
    numerical(model(device, model_name), expression)
end

function numerical(device::QuantumDevice, gate_name::Symbol)
    haskey(device.gatespecs, gate_name) ||
        error("QuantumDevice does not contain a gate named $gate_name.")
    gate = device.gatespecs[gate_name]
    numerical(model(gate.modelspec), gate)
end
