_gate_real(value) = value isa Real && !(value isa Bool) && concrete_parameter_value(value)

# Resolve aliases to real constructor fields; paths remain local to a control.
function _canonical_pulse_path(object, keys)
    isempty(keys) && return String[]
    key = first(keys)
    if object isa Pulse && length(keys)>=2 && key=="envelope" && keys[2]=="duration"
        return vcat(["duration"],String.(keys[3:end]))
    elseif object isa Segment && length(keys)>=2 && key=="shape" && keys[2]=="duration"
        return vcat(["duration"],String.(keys[3:end]))
    end
    if !hasproperty(object,Symbol(key)) && !(object isa Union{NamedTuple,Tuple,AbstractArray,AbstractDict}) && hasproperty(object,:parameters)
        return vcat(["parameters"],_canonical_pulse_path(object.parameters,keys))
    end
    child = _get(object,Symbol(key))
    return vcat([String(key)],_canonical_pulse_path(child,keys[2:end]))
end
function _gate_control_path(controls, path::AbstractString)
    keys = String.(_path_keys(path))
    length(keys) >= 2 || throw(ArgumentError("expected control_name/pulse_path: $path"))
    index = findfirst(c -> string(c.name) == first(keys),controls)
    index === nothing && throw(ArgumentError("unknown control in pulse path: $path"))
    localpath = _canonical_pulse_path(controls[index].pulse,keys[2:end])
    return index, join(vcat([string(controls[index].name)],localpath),'/'), localpath
end
function _derived_pulse_path(pulse, keys)
    isempty(keys) && return false
    key = Symbol(first(keys))
    (key === :duration && (pulse isa Union{PulseSequence,EnvelopeSequence} || pulse isa Pulse && pulse.envelope isa EnvelopeSequence || pulse isa Segment && pulse.shape isa EnvelopeSequence)) && return true
    return length(keys)>1 && _derived_pulse_path(_get(pulse,key),keys[2:end])
end
function _gate_numeric_leaf(controls,path)
    index,canonical,keys = try
        _gate_control_path(controls,path)
    catch error
        error isa KeyError || rethrow()
        throw(ArgumentError("unknown numeric pulse parameter path: $path"))
    end
    pulse = controls[index].pulse
    _derived_pulse_path(pulse,keys) && throw(ArgumentError("derived duration cannot be independently selected or linked: $path"))
    available = calibration_parameters(pulse)
    localpath = join(keys,'/')
    haskey(available,localpath) || throw(ArgumentError("not a numeric pulse parameter: $path"))
    value = getpath(pulse,localpath)
    _gate_real(value) || throw(ArgumentError("linked/selected pulse fields must be finite concrete real scalars: $path"))
    return index,canonical,keys,value
end

# Validate complete physical fields, including explicit duration ownership.
function _physical_gate(model_recipe,controls,terms,duration,duration_policy::Symbol,metadata)
    all(c -> c isa Control,controls) || throw(ArgumentError("controls must contain Controls"))
    all(t -> t isa Interaction,terms) || throw(ArgumentError("terms must contain Interactions"))
    allunique(c.name for c in controls) || throw(ArgumentError("duplicate control names"))
    allunique(c.target for c in controls) || throw(ArgumentError("duplicate control targets"))
    duration_policy in (:auto,:explicit) || throw(ArgumentError("invalid gate duration policy"))
    latest = isempty(controls) ? nothing : maximum(c.pulse.delay+c.pulse.duration for c in controls)
    finalduration = duration_policy === :auto ? latest : duration
    isnothing(finalduration) && throw(ArgumentError("duration is required without controls"))
    _require_finite_real("gate duration",finalduration)>=0 || throw(ArgumentError("duration must be nonnegative"))
    isnothing(latest) || finalduration>=latest || throw(ArgumentError("duration cuts off a control"))
    return GateRecipe(model_recipe,controls,terms,finalduration,metadata,duration_policy)
end

# Apply the entire numerical candidate in one bottom-up reconstruction per control.
function _update_gate_pulses(gate::GateRecipe,updates)
    changes = Dict{String,Any}()
    for (path,value) in updates
        _gate_real(value) || throw(ArgumentError("pulse updates must be finite concrete real scalars"))
        _,canonical,_,_ = _gate_numeric_leaf(gate.controls,path)
        haskey(changes,canonical) && throw(ArgumentError("duplicate pulse update $canonical"))
        changes[canonical] = value
    end
    controls = map(gate.controls) do control
        localchanges = [(split(path,'/')[2:end],value) for (path,value) in changes if first(split(path,'/'))==string(control.name)]
        Control(control.name,control.target,_pulse_update_tree(control.pulse,localchanges))
    end
    return _physical_gate(gate.model_recipe,controls,gate.terms,gate.duration,gate.duration_policy,gate.metadata)
end

# Structural fields retain precedence over canonical control-name pulse paths.
function getpath(gate::GateRecipe,path::AbstractString)
    keys = _path_keys(path)
    hasproperty(gate,first(keys)) && return foldl(_get,keys;init=gate)
    index,_,localkeys = _gate_control_path(gate.controls,path)
    return getpath(gate.controls[index].pulse,join(localkeys,'/'))
end
function haspath(gate::GateRecipe,path::AbstractString)
    keys = _path_keys(path)
    if !hasproperty(gate,first(keys)) && !any(c -> string(c.name)==String(first(keys)),gate.controls)
        return false
    end
    try
        getpath(gate,path)
        return true
    catch error
        error isa KeyError && return false
        rethrow()
    end
end
function setpath(gate::GateRecipe,path::AbstractString,value)
    keys = _path_keys(path)
    hasproperty(gate,first(keys)) && return _setpath(gate,keys,value)
    index,_,localkeys = _gate_control_path(gate.controls,path)
    pulse = setpath(gate.controls[index].pulse,join(localkeys,'/'),value)
    controls = [i==index ? Control(c.name,c.target,pulse) : c for (i,c) in enumerate(gate.controls)]
    return _physical_gate(gate.model_recipe,controls,gate.terms,gate.duration,gate.duration_policy,gate.metadata)
end

# Generic structural edits also revalidate physical timing and duration ownership.
function _replace_property(gate::GateRecipe,key::Symbol,value)
    hasfield(typeof(gate),key) || throw(ArgumentError("unknown GateRecipe field $key"))
    field(name) = name===key ? value : getfield(gate,name)
    policy = key===:duration ? :explicit : field(:duration_policy)
    return _physical_gate(field(:model_recipe),field(:controls),field(:terms),field(:duration),policy,field(:metadata))
end
