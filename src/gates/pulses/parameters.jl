function _path_keys(path::AbstractString)
    parts = split(path, '/')
    any(isempty, parts) && throw(ArgumentError("Paths must contain nonempty slash-separated segments"))
    return Symbol.(parts)
end

function _dict_key(object::AbstractDict, key::Symbol; insert = false)
    candidates = Any[key, String(key)]
    index = tryparse(Int, String(key))
    index === nothing || push!(candidates, index)
    matches = filter(k -> haskey(object, k), candidates)
    length(matches) > 1 && throw(ArgumentError("Ambiguous dictionary path segment: $key"))
    !isempty(matches) && return only(matches)
    insert || throw(KeyError(key))
    K = keytype(object)
    K <: AbstractString && return String(key)
    K <: Integer && index !== nothing && return index
    (K === Any || K <: Symbol) && return key
    throw(ArgumentError("Cannot create dictionary key $key of type $K"))
end

function _path_index(object, key::Symbol)
    index = tryparse(Int, String(key))
    index !== nothing && index in eachindex(object) || throw(KeyError(key))
    return index
end

_get(x::NamedTuple, key::Symbol) = haskey(x, key) ? getproperty(x, key) : throw(KeyError(key))
_get(x::AbstractDict, key::Symbol) = x[_dict_key(x, key)]
_get(x::Union{AbstractArray,Tuple}, key::Symbol) = x[_path_index(x, key)]
function _get(x, key::Symbol)
    hasproperty(x, key) && return getproperty(x, key)
    hasproperty(x, :parameters) && return _get(getproperty(x, :parameters), key)
    throw(KeyError(key))
end

_has(x::NamedTuple, key::Symbol) = hasproperty(x, key)
function _has(x, key::Symbol)
    try
        _get(x, key)
        return true
    catch error
        error isa KeyError && return false
        rethrow()
    end
end

"""
    getpath(object, path::AbstractString)

Read a nonempty slash-separated path through structs, named tuples, dictionaries,
arrays, or tuples. Array indices are one-based. A struct's parameters may omit
the `parameters/` segment when no field shadows that name. Dictionary segments
match symbol, string, or integer keys; ambiguous matches throw `ArgumentError`.
Missing keys throw `KeyError`.

```jldoctest
julia> getpath((drive=(amplitude=0.2,),), "drive/amplitude")
0.2
```
"""
getpath(object, path::AbstractString) = foldl(_get, _path_keys(path); init = object)
"""
    haspath(object, path::AbstractString)

Return whether a path resolves using [`getpath`](@ref) rules. Missing keys return
`false`; malformed or ambiguous paths still throw rather than hiding the error.
"""
function haspath(object, path::AbstractString)
    for key in _path_keys(path)
        _has(object, key) || return false
        object = _get(object, key)
    end
    return true
end
# Rebuild through the outer constructor, allowing parametric field types to change.
# Domain types can specialize this hook next to their own definitions.
function _replace_property(object, key::Symbol, value)
    hasfield(typeof(object), key) || throw(ArgumentError("Property $key is not a writable field"))
    fields = map(fieldnames(typeof(object))) do name
        name === key ? value : getfield(object, name)
    end
    constructor = Base.typename(typeof(object)).wrapper
    return constructor(fields...)
end

function _replace_property(object::NamedTuple, key::Symbol, value)
    haskey(object, key) || throw(KeyError(key))
    return merge(object, NamedTuple{(key,)}((value,)))
end
function _replace_property(object::Tuple, key::Symbol, value)
    index = _path_index(object, key)
    return ntuple(i -> i == index ? value : object[i], length(object))
end
function _replace_property(object::AbstractArray, key::Symbol, value)
    index = _path_index(object, key)
    T = promote_type(eltype(object), typeof(value))
    result = similar(object, T)
    copyto!(result, object)
    result[index] = value
    return result
end
function _replace_property(object::AbstractDict, key::Symbol, value)
    actual_key = _dict_key(object, key; insert = true)
    result = copy(object)
    result[actual_key] = value
    return result
end

function _setpath(object, keys, value)
    key = first(keys)
    # Match getpath's shorthand lookup without treating computed properties as fields.
    if !(object isa Union{NamedTuple,Tuple,AbstractArray,AbstractDict}) && !hasproperty(object, key)
        hasproperty(object, :parameters) || throw(KeyError(key))
        return _replace_property(object, :parameters, _setpath(object.parameters, keys, value))
    end
    replacement = length(keys) == 1 ? value : _setpath(_get(object, key), keys[2:end], value)
    return _replace_property(object, key, replacement)
end

"""
    setpath(object, path::AbstractString, value)

Return an updated object by copying containers along the path and reconstructing
structs. Unchanged branches may be shared; the original is not modified.
```jldoctest
julia> setpath((amplitude=0.2,), "amplitude", 0.4)
(amplitude = 0.4,)
```
"""
setpath(object, path::AbstractString, value) = _setpath(object, _path_keys(path), value)

const _ShapePayload=Union{Constant,Gaussian,GaussianZero,SineSquared,Cosine,Blackman,Sech,Bump,DRAG,Slepian,GaussianSquare,ErfSquare,RampedFlattop,InterpolatedEnvelope,CustomShape,Region}
Envelope(s::_ShapePayload;duration,rise=nothing,fall=nothing)=_wrap_shape(s,duration,rise,fall)
const _PulseObject = Union{AbstractPulse,AbstractEnvelope,ConstantFrequency,Segment,_ShapePayload}
function _pulse_fields(x)
    x isa Envelope && return (:duration,keys(_physical_parameters(x.shape,x.duration))...)
    x isa LegacyEnvelope && return (:duration,keys(x.fields)...)
    x isa Pulse && return (:envelope,:carrier,:amplitude,:offset,:duration,:delay)
    fields=fieldnames(typeof(x))
    return filter(k -> !(k in (:interpolant, :weights)) && !(x isa PulseSequence && k===:duration), fields)
end
function _parameter_leaves!(result, x, path;skip_duration=false)
    x isa _ShapePayload && return result
    if x isa _PulseObject || x isa NamedTuple
        for k in (x isa NamedTuple ? keys(x) : _pulse_fields(x))
            skip_duration && k===:duration && continue
            owned_child=(x isa Pulse && k===:envelope) || (x isa Segment && k===:shape)
            _parameter_leaves!(result,getproperty(x,k),isempty(path) ? string(k) : "$path/$k";skip_duration=owned_child)
        end
    elseif x isa Union{Tuple,AbstractVector}
        for i in eachindex(x)
            _parameter_leaves!(result, x[i], "$path/$i")
        end
    else
        result[path]=(path, x)
    end
    return result
end
"""All explicit pulse fields as slash paths and `(path,value)` pairs, including structural choices."""
parameters(x::_ShapePayload)=throw(ArgumentError("Normalized payload fields are internal; discover physical parameters on its Envelope"))
parameters(x::_PulseObject) = _parameter_leaves!(Dict{String,Tuple}(), x, "")
"""Discover numeric candidates and symbolic offset references; excludes Boolean and sample-count configuration. Symbolic references are not numeric optimizer variables."""
function calibration_parameters(x::_PulseObject)
    return filter(parameters(x)) do entry
        path, (_, value) = entry
        parts = split(path, '/')
        parent = length(parts) == 1 ? x : getpath(x, join(parts[1:(end-1)], '/'))
        value isa Number && !(value isa Bool) && !(parent isa Envelope && parent.shape isa Slepian && last(parts) == "samples")
    end
end
_replace_property(x::_PulseObject,k::Symbol,v)=_rebuild_pulse_node(x,NamedTuple{(k,)}((v,)))
setpath(x::_PulseObject,path::AbstractString,value)=_pulse_update_tree(x,[(_path_keys(path),value)])

# Reconstruct a node from all final fields, rather than validating intermediate edits.
function _rebuild_pulse_node(object, replacements::NamedTuple)
    value(key) = haskey(replacements,key) ? getproperty(replacements,key) : getproperty(object,key)
    if object isa Envelope
        d=value(:duration)
        allowed=(:duration,:shape,:rise,:fall,keys(_physical_parameters(object.shape,d))...)
        all(k->k in allowed,keys(replacements)) || throw(ArgumentError("Unknown or structural envelope parameter"))
        shape=get(replacements,:shape,object.shape)
        edits=Base.structdiff(replacements,(duration=nothing,shape=nothing,rise=nothing,fall=nothing))
        if !isempty(edits)
            shape=_shape_from_physical(shape,merge(_physical_parameters(shape,d),edits),d)
        end
        return _raw_envelope(shape,d,get(replacements,:rise,object.rise),get(replacements,:fall,object.fall))
    elseif object isa LegacyEnvelope
        d=get(replacements,:duration,object.duration)
        edits=Base.structdiff(replacements,(duration=nothing,))
        isempty(edits) && d!==nothing && return _stretch_envelope(object,d)
        return LegacyEnvelope(object.implementation,merge(object.fields,edits),d)
    elseif object isa Pulse
        e = value(:envelope)
        d = haskey(replacements,:duration) ? replacements.duration : _envelope_duration(e)
        return Pulse(e,value(:carrier),value(:amplitude),value(:offset),d,value(:delay))
    elseif object isa PulseSequence
        all(k->k in (:pulses,:offset,:delay),keys(replacements)) || throw(ArgumentError("PulseSequence duration is derived"))
        return PulseSequence(value(:pulses)...;offset=value(:offset),delay=value(:delay))
    elseif object isa EnvelopeSequence
        all(k->k===:segments,keys(replacements)) || throw(ArgumentError("EnvelopeSequence duration is derived; edit segment durations"))
        return EnvelopeSequence(value(:segments)...)
    elseif object isa Segment
        e=value(:shape)
        d=e isa EnvelopeSequence && !haskey(replacements,:duration) ? _envelope_duration(e) : value(:duration)
        return Segment(e;duration=d,scale=value(:scale))
    elseif object isa NamedTuple
        return merge(object,replacements)
    end
    constructor = Base.typename(typeof(object)).wrapper
    return constructor((value(k) for k in fieldnames(typeof(object)))...)
end

function _pulse_update_tree(object, changes)
    isempty(changes) && return object
    roots = filter(c -> isempty(first(c)),changes)
    if !isempty(roots)
        length(changes) == 1 || throw(ArgumentError("overlapping pulse updates"))
        return last(only(roots))
    end
    object isa _ShapePayload && throw(ArgumentError("Normalized shape payloads are read-only; edit physical fields on Envelope"))
    grouped = Dict{Symbol,Vector{Any}}()
    for (path,value) in changes
        push!(get!(grouped,Symbol(first(path)),Any[]),(path[2:end],value))
    end
    if object isa Envelope && haskey(grouped,:duration)
        entries=grouped[:duration]
        length(entries)==1 && isempty(first(only(entries))) || throw(ArgumentError("Invalid duration edit"))
        object=_stretch_envelope(object,last(only(entries)))
    end
    if object isa Union{Tuple,AbstractArray}
        entries = Dict(_path_index(object,k)=>_pulse_update_tree(_get(object,k),entries)
            for (k,entries) in grouped)
        object isa Tuple && return ntuple(i -> get(entries,i,object[i]),length(object))
        T = promote_type(eltype(object),(typeof(v) for v in values(entries))...)
        result = similar(object,T)
        copyto!(result,object)
        for (i,v) in entries
            result[i] = v
        end
        return result
    elseif object isa AbstractDict
        result = copy(object)
        for (key,entries) in grouped
            actual = _dict_key(object,key)
            result[actual] = _pulse_update_tree(object[actual],entries)
        end
        return result
    end
    if object isa Segment && haskey(grouped,:shape)
        durations=filter(c->first(c)==[:duration] || first(c)==["duration"],grouped[:shape])
        if !isempty(durations)
            haskey(grouped,:duration) && throw(ArgumentError("Duplicate duration aliases"))
            grouped[:duration]=Any[([],last(only(durations)))]
            filter!(c->!(first(c)==[:duration] || first(c)==["duration"]),grouped[:shape])
        end
    end
    if object isa Union{Pulse,Segment} && haskey(grouped,:duration)
        childkey=object isa Pulse ? :envelope : :shape
        child=getproperty(object,childkey)
        if child isa Envelope
            dchanges=grouped[:duration]
            length(dchanges)==1 && isempty(first(only(dchanges))) || throw(ArgumentError("Invalid duration edit"))
            entries=get!(grouped,childkey,Any[])
            any(c->first(c)==[:duration] || first(c)==["duration"],entries) && throw(ArgumentError("Duplicate duration aliases"))
            push!(entries,([:duration],last(only(dchanges))))
        end
    end
    # Normalize custom Envelope shorthand before grouping/reconstructing.
    direct = Pair{Symbol,Any}[]
    parameter_changes = Any[]
    for (key,entries) in grouped
        if !hasproperty(object,key) && hasproperty(object,:parameters)
            append!(parameter_changes,[(vcat([String(key)],path),v) for (path,v) in entries])
        else
            push!(direct,key=>_pulse_update_tree(_get(object,key),entries))
        end
    end
    if !isempty(parameter_changes)
        any(p -> first(p) === :parameters,direct) && throw(ArgumentError("ambiguous parameter aliases"))
        push!(direct,:parameters=>_pulse_update_tree(object.parameters,parameter_changes))
    end
    return _rebuild_pulse_node(object,(;direct...))
end
