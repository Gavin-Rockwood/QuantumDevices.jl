# Known complex envelope payloads need a carrier's real laboratory-frame mixer.
# Custom callables remain responsible for real finite values over their domain.
function _validate_voltage_shape(shape,tree)
    shape isa Union{DRAG,_LegacyPhysical.DRAG} && throw(ArgumentError("complex DRAG voltage envelope requires a real carrier mixer"))
    if shape isa Union{InterpolatedEnvelope,_LegacyPhysical.InterpolatedEnvelope}
        all(_map_real,tree.values) || throw(ArgumentError("unmodulated voltage envelope samples must be finite real volts"))
    elseif shape isa Region
        _validate_voltage_shape(shape.shape,tree)
    elseif shape isa RampedFlattop
        for key in (:ramp_up,:ramp_down)
            ramp = getproperty(shape,key)
            ramp isa NamedTuple && _validate_voltage_shape(ramp.shape,getproperty(tree,key))
        end
    end
    return nothing
end
function _validate_voltage_envelope(envelope,tree)
    if envelope isa Envelope
        _validate_voltage_shape(envelope.shape,tree)
    elseif envelope isa EnvelopeSequence
        for (segment,values) in zip(envelope.segments,tree.segments)
            _map_real(values.scale) || throw(ArgumentError("unmodulated voltage segment scales must be finite real values"))
            _validate_voltage_envelope(segment.shape,values.shape)
        end
    elseif envelope isa LegacyEnvelope
        _validate_voltage_shape(_legacy_value(envelope),tree)
    elseif envelope isa PhysicalEnvelope
        _validate_voltage_envelope(envelope.source,tree)
    end
    return nothing
end
# Mapped voltage offsets and amplitudes are explicit numeric values at every level.
function _validate_voltage_tree(recipe,tree)
    _map_real(tree.offset) || throw(ArgumentError("mapped voltage offsets must be explicit finite real volts"))
    if recipe isa PulseSequenceRecipe
        foreach(pair -> _validate_voltage_tree(pair...),zip(recipe.pulses,tree.pulses))
    else
        _map_real(tree.amplitude) || throw(ArgumentError("mapped voltage amplitudes must be finite real volts"))
        recipe.carrier === nothing && _validate_voltage_envelope(recipe.envelope,tree.envelope)
    end
    return nothing
end
_validate_voltage_pulse(recipe) = _validate_voltage_tree(recipe,_parameter_tree(recipe))
# Resolve routing once at a public setup boundary; no route lookup occurs in callbacks.
function _prepare_control_routes(gate,values,required)
    direct = Control[]
    groups = Pair{AbstractControlMap,Tuple}[]
    selected = Symbol[]
    targets = Set{Num}()
    maps = gate.model_recipe isa Recipe ? gate.model_recipe.control_maps : throw(ArgumentError("capture the named gate endpoint before routing"))
    for c in gate.controls
        if c.target isa Num
            c.target in targets && throw(ArgumentError("duplicate physical control target"))
            push!(targets,c.target)
            push!(direct,c)
        else
            address = c.target
            haskey(maps,address.mapname) || throw(ArgumentError("unknown captured control map $(address.mapname)"))
            address.inputname in maps[address.mapname].inputs || throw(ArgumentError("unknown map input $(address.inputname)"))
            _validate_voltage_pulse(c.pulse.recipe)
            address.mapname in selected || push!(selected,address.mapname)
        end
    end
    for name in selected
        map = maps[name]
        controls = Base.map(map.inputs) do input
            matches = [c for c in gate.controls if c.target isa ControlInput && isequal(c.target,ControlInput(name,input))]
            length(matches)==1 || throw(ArgumentError("map $name requires exactly one control for input $input"))
            only(matches)
        end
        isempty(intersect(targets,Set(map.outputs))) || throw(ArgumentError("direct/mapped or mapped/mapped physical output collision"))
        union!(targets,map.outputs)
        push!(groups,map=>controls)
    end
    issubset(targets,required) || throw(ArgumentError("control target must appear in endpoint Hamiltonian or gate terms"))
    return (;direct,groups,targets)
end
# Candidate rebuilding and waveform compilation stay outside the solver RHS.
function _candidate_control_callbacks(gate,reference,required)
    bindings = _compile_gate_bindings(gate,reference,();required)
    return Dict(variable => let f=callback; t -> f(nothing,t); end for (variable,callback) in bindings)
end
"""Sample map domains outside evolution; finite samples do not guarantee continuous coverage."""
function control_map_diagnostics(gate::GateRecipe,model::NumericalModel;times,params=nothing)
    samples = collect(times)
    all(_map_real,samples) || throw(ArgumentError("diagnostic times must be finite real values"))
    _, values = _reuse_gate_data(model,gate)
    required = union(refs(model.hamiltonian,:param),Set{Num}(v for term in gate.terms for v in refs(term.expression,:param)))
    routes = _prepare_control_routes(gate,values,required)
    paths = isnothing(params) ? () : Tuple(keys(params))
    _prepare_gate_runtime(gate,values,paths,params)
    selected = _runtime_gate_paths(gate,paths)
    reports = map(routes.groups) do (map,controls)
        waveforms = Tuple(_compiled_control_waveform(c,values,selected) for c in controls)
        inputs = [Tuple(f(params,t) for f in waveforms) for t in samples]
        checks = Base.map(inputs) do input
            try
                (;output=check_control_map(map,input),failure=nothing)
            catch error
                error isa Union{ArgumentError,DomainError} || rethrow()
                (;output=nothing,failure=sprint(showerror,error))
            end
        end
        (;name=map.name,inputs,outputs=[c.output for c in checks],failures=[c.failure for c in checks])
    end
    return (;times=samples,maps=reports,sampled_only=true)
end
