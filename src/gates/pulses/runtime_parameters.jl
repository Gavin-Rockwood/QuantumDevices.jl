# Public overrides and solve preparation resolve and validate final fields.
# Public keyword calls resolve explicit physical edits at the call boundary.
function _compile_pulse(recipe,defaults)
    direct=_compile_waveform(recipe)
    allpaths=Tuple(path for path in keys(_numeric_parameter_dict(defaults)) if !_structural_pulse_parameter(recipe,split(path,'/')) && !_rebuild_only_parameter(recipe,split(path,'/')))
    physical=_compile_waveform(recipe,allpaths)
    return function(t;kwargs...)
        _require_finite_real("time",t)
        if isempty(kwargs)
            _require_resolved_pulse(recipe)
            return direct(nothing,t,t)
        end
        flat=Dict{String,Any}()
        nested=Pair{Symbol,Any}[]
        for (key,value) in pairs(kwargs)
            if occursin('/',String(key))
                haskey(_numeric_parameter_dict(defaults),String(key)) || throw(ArgumentError("unknown pulse parameter $key"))
                flat[String(key)]=value
            else
                push!(nested,key=>value)
            end
        end
        changes=_keyword_changes(defaults,(;nested...))
        supplied=Dict{String,Any}(join(string.(path),'/')=>value for (path,value) in changes)
        isempty(intersect(keys(supplied),keys(flat))) || throw(ArgumentError("duplicate physical pulse parameter through flat/nested aliases"))
        merge!(supplied,flat)
        original=_numeric_parameter_dict(defaults)
        for (path,value) in supplied
            haskey(original,path) && !_structural_pulse_parameter(recipe,split(path,'/')) || throw(ArgumentError("unknown or structural pulse parameter $path"))
            _derived_pulse_path(recipe,split(path,'/')) && throw(ArgumentError("derived pulse duration is read-only: $path"))
            if _rebuild_only_parameter(recipe,split(path,'/')) && !isequal(value,original[path])
                throw(ArgumentError("$path requires rebuilding; use setpath(pulse,path,value)"))
            end
        end
        resolved=_resolve_physical_parameters(recipe,defaults,supplied)
        return physical(_concrete_runtime_dict(_numeric_parameter_dict(resolved)),t,t)
    end
end
function _keyword_changes(defaults,overrides,path=Symbol[])
    changes=Any[]
    for (key,value) in pairs(overrides)
        if defaults isa NamedTuple
            haskey(defaults,key) || throw(ArgumentError("unknown pulse parameter $key"))
        elseif defaults isa Tuple
            key=Symbol(key)
        end
        original=_get(defaults,key)
        next=vcat(path,[key])
        if value isa NamedTuple && original isa NamedTuple
            append!(changes,_keyword_changes(original,value,next))
        elseif value isa Tuple && original isa Tuple
            length(value)==length(original) || throw(ArgumentError("pulse tuple size is structural"))
            append!(changes,_keyword_changes(original,(; (Symbol(i)=>v for (i,v) in enumerate(value))...),next))
        else
            push!(changes,(next,value))
        end
    end
    return changes
end
function _numeric_parameter_dict(tree,base="",out=Dict{String,Any}())
    if tree isa NamedTuple
        for (k,v) in pairs(tree); _numeric_parameter_dict(v,_path_join(base,k),out); end
    elseif tree isa Tuple
        for (i,v) in enumerate(tree); _numeric_parameter_dict(v,_path_join(base,i),out); end
    elseif tree isa Number
        out[base]=tree
    end
    return out
end
function _resolve_tree(defaults,updates,base="";ratio=1)
    if defaults isa NamedTuple
        return NamedTuple{keys(defaults)}(map(keys(defaults)) do k
            value=getproperty(defaults,k)
            factor=k in (:sigma,:beta,:width,:center,:ramp_time,:rise_time,:fall_time,:times) ? ratio : 1
            _resolve_tree(value,updates,_path_join(base,k);ratio=factor)
        end)
    elseif defaults isa Tuple
        return map(Tuple(enumerate(defaults))) do (i,v)
            _resolve_tree(v,updates,_path_join(base,i);ratio)
        end
    end
    return get(updates,base,defaults isa Number ? defaults*ratio : defaults)
end
function _resolve_physical_parameters(r::PulseRecipe,defaults,updates;base="")
    duration=get(updates,_path_join(base,:duration),defaults.duration)
    envelope=if r.envelope isa EnvelopeSequence
        _resolve_envelope_sequence(r.envelope,defaults.envelope,updates,_path_join(base,:envelope))
    else
        r.envelope isa Envelope ? _resolve_shape_parameters(r.envelope.shape,defaults.envelope,updates,_path_join(base,:envelope),duration,defaults.duration) : _resolve_tree(defaults.envelope,updates,_path_join(base,:envelope))
    end
    if r.envelope isa EnvelopeSequence
        duration=foldl(_checked_time_add,(s.duration for s in envelope.segments))
    end
    result=(amplitude=get(updates,_path_join(base,:amplitude),defaults.amplitude),
        offset=get(updates,_path_join(base,:offset),defaults.offset),duration=duration,
        delay=get(updates,_path_join(base,:delay),defaults.delay),
        carrier=_resolve_tree(defaults.carrier,updates,_path_join(base,:carrier)),envelope=envelope)
    _validate_resolved_pulse(r,result)
    return result
end
function _resolve_physical_parameters(r::PulseSequenceRecipe,defaults,updates;base="")
    children=map(Tuple(eachindex(r.pulses))) do i
        _resolve_physical_parameters(r.pulses[i],defaults.pulses[i],updates;base=_path_join(base,"pulses/$i"))
    end
    result=(pulses=children,offset=get(updates,_path_join(base,:offset),defaults.offset),delay=get(updates,_path_join(base,:delay),defaults.delay))
    _require_finite("offset",result.offset)
    _require_finite_real("delay",result.delay)>=0 || throw(ArgumentError("delay must be nonnegative"))
    return result
end
function _resolve_envelope_sequence(e,defaults,updates,base)
    segments=map(Tuple(eachindex(e.segments))) do i
        s=e.segments[i]; old=defaults.segments[i]; prefix=_path_join(base,"segments/$i")
        d=get(updates,_path_join(prefix,:duration),old.duration)
        shape=if s.shape isa EnvelopeSequence
            _resolve_envelope_sequence(s.shape,old.shape,updates,_path_join(prefix,:shape))
        else
            s.shape isa Envelope ? _resolve_shape_parameters(s.shape.shape,old.shape,updates,_path_join(prefix,:shape),d,old.duration) : _resolve_tree(old.shape,updates,_path_join(prefix,:shape))
        end
        s.shape isa EnvelopeSequence && (d=foldl(_checked_time_add,(child.duration for child in shape.segments)))
        _require_positive_width("segment duration",d)
        scale=get(updates,_path_join(prefix,:scale),old.scale); _require_finite("segment scale",scale)
        _validate_resolved_envelope(s.shape,shape,d)
        (shape=shape,duration=d,scale=scale)
    end
    return (segments=segments,)
end
function _validate_resolved_pulse(r::PulseRecipe,p)
    _require_positive_width("duration",p.duration)
    _require_finite_real("delay",p.delay)>=0 || throw(ArgumentError("delay must be nonnegative"))
    _checked_time_add(p.delay,p.duration)
    _require_finite("amplitude",p.amplitude); _require_finite("offset",p.offset)
    if p.carrier!==nothing
        _require_finite_real("amplitude",p.amplitude); _require_finite_real("offset",p.offset)
        _require_finite_real("frequency",p.carrier.frequency); _require_finite_real("phase",p.carrier.phase)
    end
    _validate_resolved_envelope(r.envelope,p.envelope,p.duration)
    return nothing
end
_validate_resolved_envelope(e::Envelope,p,d) = _validate_resolved_shape(e.shape,p,d)
function _validate_resolved_shape(s,p,d)
    if s isa Region
        return _validate_resolved_shape(s.shape,p,d/(s.interval[2]-s.interval[1]))
    end
    for k in (:sigma,:width)
        haskey(p,k) && _require_positive_width(string(k),getproperty(p,k))
    end
    haskey(p,:beta) && _require_finite_real("beta",p.beta)
    if haskey(p,:center)
        _require_finite_real("center",p.center)
        0<=p.center<=d || throw(ArgumentError("center must lie inside duration"))
    end
    haskey(p,:k) && _require_positive_width("k",p.k)
    if s isa Union{GaussianSquare,ErfSquare}
        _require_finite_real("ramp_time",p.ramp_time)
        0<=p.ramp_time<=d/2 || throw(ArgumentError("both ramps must fit duration"))
    elseif s isa RampedFlattop
        _require_finite_real("rise_time",p.rise_time); _require_finite_real("fall_time",p.fall_time)
        0<=p.rise_time && 0<=p.fall_time && p.rise_time+p.fall_time<=d || throw(ArgumentError("ramps must fit duration"))
        s.ramp_up isa NamedTuple && p.rise_time>0 && _validate_resolved_shape(s.ramp_up.shape,p.ramp_up,p.rise_time)
        s.ramp_down isa NamedTuple && p.fall_time>0 && _validate_resolved_shape(s.ramp_down.shape,p.ramp_down,p.fall_time)
    elseif s isa InterpolatedEnvelope
        foreach(t -> _require_finite_real("sample time",t),p.times)
        first(p.times)==0 && last(p.times)==d && all(i -> p.times[i]<p.times[i+1],1:length(p.times)-1) || throw(ArgumentError("sample nodes must increase from zero to duration"))
        all(isfinite,p.times) && all(isfinite,p.values) || throw(ArgumentError("nonfinite samples"))
    end
    return nothing
end
_validate_resolved_envelope(e::EnvelopeSequence,p,d)=nothing
function _resolve_reference_tree(r::PulseRecipe,tree,reference)
    return merge(tree,(offset=_reference_offset(tree.offset,reference),))
end
function _resolve_reference_tree(r::PulseSequenceRecipe,tree,reference)
    return merge(tree,(offset=_reference_offset(tree.offset,reference),pulses=map((child,v)->_resolve_reference_tree(child,v,reference),r.pulses,tree.pulses)))
end
_base_shape(s::Region)=_base_shape(s.shape)
_base_shape(s)=s
function _rebuild_only_parameter(object,keys)
    object isa Pulse && return _rebuild_only_parameter(object.recipe,keys)
    isempty(keys) && return false
    key=Symbol(first(keys))
    if object isa Envelope
        shape=_base_shape(object.shape)
        shape isa Slepian && key===:time_bandwidth && return true
        shape isa CustomShape && key===:parameters && return true
    elseif object isa AbstractEnvelope && !(object isa Union{Envelope,EnvelopeSequence})
        return true
    end
    return length(keys)>1 && _rebuild_only_parameter(_get(object,key),keys[2:end])
end
_validate_resolved_envelope(e::LegacyEnvelope,p,d)=nothing
function _resolve_shape_parameters(s,defaults,updates,base,d,old_d)
    s isa Region && return _resolve_shape_parameters(s.shape,defaults,updates,base,d/(s.interval[2]-s.interval[1]),old_d/(s.interval[2]-s.interval[1]))
    resolved=_resolve_tree(defaults,updates,base;ratio=d/old_d)
    if s isa RampedFlattop
        if s.ramp_up isa NamedTuple && defaults.ramp_up isa NamedTuple
            source_defaults=iszero(s.rise_time) ? _physical_tree(_physical_parameters(s.ramp_up.shape,one(old_d))) : defaults.ramp_up
            source_duration=iszero(s.rise_time) ? one(old_d) : s.rise_time*old_d
            up=_resolve_shape_parameters(s.ramp_up.shape,source_defaults,updates,_path_join(base,:ramp_up),resolved.rise_time,source_duration)
            resolved=merge(resolved,(ramp_up=up,))
        end
        if s.ramp_down isa NamedTuple && defaults.ramp_down isa NamedTuple
            source_defaults=iszero(s.fall_time) ? _physical_tree(_physical_parameters(s.ramp_down.shape,one(old_d))) : defaults.ramp_down
            source_duration=iszero(s.fall_time) ? one(old_d) : s.fall_time*old_d
            down=_resolve_shape_parameters(s.ramp_down.shape,source_defaults,updates,_path_join(base,:ramp_down),resolved.fall_time,source_duration)
            resolved=merge(resolved,(ramp_down=down,))
        end
    end
    return resolved
end

function _concrete_runtime_dict(entries)
    isempty(entries) && return Dict{String,Float64}()
    types=Tuple(unique(typeof(value) for value in values(entries)))
    T=all(value -> value isa Real,values(entries)) || all(value -> value isa Complex,values(entries)) ? promote_type(types...) : Union{types...}
    return Dict{String,T}(String(key)=>value for (key,value) in pairs(entries))
end
_validate_resolved_envelope(e::AbstractEnvelope,p,d)=validate_envelope(e,d)

function _structural_pulse_parameter(object,keys)
    object isa Pulse && return _structural_pulse_parameter(object.recipe,keys)
    isempty(keys) && return false
    key=Symbol(first(keys))
    object isa Envelope && _base_shape(object.shape) isa Slepian && key===:samples && return true
    return length(keys)>1 && _structural_pulse_parameter(_get(object,key),keys[2:end])
end
