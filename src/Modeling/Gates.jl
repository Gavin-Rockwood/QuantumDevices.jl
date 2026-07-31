"""A slow envelope multiplied by a fixed laboratory-frame carrier."""
struct CarrierControl{E}
    envelope::E
    frequency::Float64
    phase::Float64
end

function CarrierControl(envelope; frequency::Real, phase::Real = 0.0)
    isfinite(frequency) && frequency > 0 ||
        error("CarrierControl frequency must be finite and positive.")
    isfinite(phase) || error("CarrierControl phase must be finite.")
    CarrierControl(envelope, Float64(frequency), Float64(phase))
end

_carrier_modulation(control::CarrierControl) =
    time -> cos(2π * control.frequency * time + control.phase)

"""A finite-duration binding of one model's exposed Hamiltonian controls."""
struct GateSpec{M<:ModelSpec,P<:NamedTuple}
    name::Symbol
    modelspec::M
    duration::Float64
    controls::Dict{ParamPath,Any}
    parameters::P
    control_recipes::Union{Nothing,Dict{ParamPath,Any}}
end

function GateSpec(
    name::Symbol,
    modelspec::ModelSpec;
    duration::Real,
    parameters = nothing,
    controls = Dict{ParamPath,Any}(),
)
    duration > 0 || error("GateSpec duration must be positive.")
    controls isa AbstractDict ||
        error("GateSpec controls must be an AbstractDict with ParamPath keys.")

    prepared_parameters = _gate_parameters(parameters)
    context = merge(prepared_parameters, (; duration = Float64(duration)))
    recipes = parameters === nothing ? nothing : Dict{ParamPath,Any}()
    active = _hamiltonian_parameters(modelspec)
    registry = Dict{ParamPath,Any}()
    for (path, value) in controls
        path isa ParamPath ||
            error("GateSpec control keys must be ParamPath objects; got $(typeof(path)).")
        haskey(active, path) || begin
            if haskey(modelspec.parameters, path) && modelspec.parameters[path].fixed
                error("GateSpec control $(join(path.parts, "/")) is fixed.")
            end
            error(
                "GateSpec control $(join(path.parts, "/")) is not exposed by its ModelSpec.",
            )
        end
        resolved = _resolve_gate_control(path, value, context, parameters !== nothing)
        checked = _parameter_value(active[path], resolved)
        _check_control_endpoint(modelspec, path, checked, 0.0)
        _check_control_endpoint(modelspec, path, checked, Float64(duration))
        registry[path] = checked
        recipes === nothing || (recipes[path] = value)
    end
    GateSpec(name, modelspec, Float64(duration), registry, prepared_parameters, recipes)
end

function _gate_parameters(values)
    values === nothing && return NamedTuple()
    values isa NamedTuple || error("GateSpec parameters must be a NamedTuple.")
    hasproperty(values, :duration) &&
        error("GateSpec parameters must not contain the reserved name :duration.")
    for (name, value) in pairs(values)
        value isa Real || error("GateSpec parameter $name must be a real scalar.")
        isfinite(value) || error("GateSpec parameter $name must be finite.")
    end
    values
end

function _resolve_gate_control(path, value, context, parameterized)
    if value isa CarrierControl
        value.envelope isa Number && error(
            "CarrierControl envelope for $(join(path.parts, "/")) must accept " *
            (parameterized ? "(t) or (parameters, t)." : "(t)."),
        )
        envelope = _resolve_gate_callable(
            path,
            value.envelope,
            context,
            parameterized;
            label = "CarrierControl envelope",
        )
        modulation = _carrier_modulation(value)
        return time -> envelope(time) * modulation(time)
    end
    value isa Number && return value
    _resolve_gate_callable(path, value, context, parameterized; label = "GateSpec control")
end

function _resolve_gate_callable(path, value, context, parameterized; label)
    accepts_time = applicable(value, 0.0)
    accepts_parameters = parameterized && applicable(value, context, 0.0)
    accepts_time &&
        accepts_parameters &&
        error(
            "$label $(join(path.parts, "/")) has ambiguous callable arity; " *
            "define either (t) or (parameters, t), not both.",
        )
    accepts_parameters && return time -> value(context, time)
    accepts_time && return value
    expected = parameterized ? "(t) or (parameters, t)" : "(t)"
    error("$label $(join(path.parts, "/")) must accept $expected.")
end

function _carrier_envelope(gate::GateSpec, path::ParamPath)
    gate.control_recipes === nothing && return nothing
    recipe = get(gate.control_recipes, path, nothing)
    recipe isa CarrierControl || return nothing
    context = merge(gate.parameters, (; duration = gate.duration))
    envelope = _resolve_gate_callable(
        path,
        recipe.envelope,
        context,
        true;
        label = "CarrierControl envelope",
    )
    (; envelope, modulation = _carrier_modulation(recipe), carrier = recipe)
end

"""Return an independently reconstructed gate with updated named recipe values."""
function with_parameters(gate::GateSpec; kwargs...)
    gate.control_recipes === nothing &&
        error("GateSpec $(gate.name) does not contain a parameterized control recipe.")
    updates = (; kwargs...)
    allowed = (propertynames(gate.parameters)..., :duration)
    unknown = filter(name -> name ∉ allowed, propertynames(updates))
    if !isempty(unknown)
        suffix = length(unknown) == 1 ? "" : "s"
        error("Unknown GateSpec parameter$suffix: " * join(string.(unknown), ", "))
    end
    duration = hasproperty(updates, :duration) ? updates.duration : gate.duration
    parameter_updates = NamedTuple{Tuple(filter(!=(:duration), propertynames(updates)))}(
        Tuple(updates[name] for name in filter(!=(:duration), propertynames(updates))),
    )
    parameters = merge(gate.parameters, parameter_updates)
    GateSpec(
        gate.name,
        gate.modelspec;
        duration,
        parameters,
        controls = deepcopy(gate.control_recipes),
    )
end

function _check_control_endpoint(spec, path, value, time)
    actual = _at(value, time)
    expected = _at(_model_value(spec, path), time)
    isapprox(actual, expected; atol = 1e-10, rtol = 1e-10) || error(
        "GateSpec control $(join(path.parts, "/")) must return to its " *
        "model default at t=$time; got $actual instead of $expected.",
    )
end

"""Build the complete static or time-dependent Hamiltonian for a gate."""
function numerical(model::QuantumDeviceModel, gate::GateSpec)
    model.spec === gate.modelspec ||
        error("GateSpec $(gate.name) belongs to a different ModelSpec.")
    _bound_model_hamiltonian(model, gate.controls)
end

numerical(gate::GateSpec) = numerical(model(gate.modelspec), gate)

function _snapshot_gate(gate::GateSpec)
    snapshot = deepcopy(gate.modelspec)
    if gate.control_recipes === nothing
        return GateSpec(
            gate.name,
            snapshot;
            duration = gate.duration,
            controls = deepcopy(gate.controls),
        )
    end
    GateSpec(
        gate.name,
        snapshot;
        duration = gate.duration,
        parameters = deepcopy(gate.parameters),
        controls = deepcopy(gate.control_recipes),
    )
end
