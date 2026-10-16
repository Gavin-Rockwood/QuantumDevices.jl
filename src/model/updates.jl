# Canonical update trees let multiple leaves reach a domain constructor together.
# They are prepared once for tomography; ordinary setpath uses the same domain
# reconstruction hooks for its single update.
mutable struct _ModelUpdateNode
    index::Int
    children::Dict{Symbol,_ModelUpdateNode}
end
_ModelUpdateNode() = _ModelUpdateNode(0, Dict{Symbol,_ModelUpdateNode}())

function _check_update_key(object, key)
    hasfield(typeof(object), key) || throw(ArgumentError("Property $key is derived or unknown"))
end
_check_update_key(::Union{NamedTuple,AbstractDict,AbstractArray,Tuple}, key) = nothing
function _check_update_key(::DeviceModel, key)
    key in _MODEL_INPUT_FIELDS || throw(ArgumentError("DeviceModel.$key is derived or unknown"))
end
function _check_update_key(::Component, key)
    key in (:name, :parameters, :dimension) || throw(ArgumentError("Component.$key is derived or unknown"))
end

function _canonical_update_path(object, keys)
    key = first(keys)
    if !(object isa Union{NamedTuple,AbstractDict,AbstractArray,Tuple}) && !hasproperty(object, key)
        hasproperty(object, :parameters) || throw(KeyError(key))
        return _canonical_update_path(object, vcat([:parameters], keys))
    end
    _check_update_key(object, key)
    if object isa Component && key === :parameters && length(keys) > 1 &&
            object.type == "TunableTransmon" && keys[2] in (:EJmax, :d)
        throw(ArgumentError("$(keys[2]) is derived from physical parameters"))
    end
    if object isa AbstractDict
        key = Symbol(string(_dict_key(object, key)))
    elseif object isa Union{AbstractArray,Tuple}
        key = Symbol(string(_path_index(object, key)))
    end
    child = _get(object, key) # Also checks dictionary ambiguity and array bounds.
    return length(keys) == 1 ? [key] : vcat([key], _canonical_update_path(child, keys[2:end]))
end

function _prepare_model_updates(model::DeviceModel, paths)
    root = _ModelUpdateNode()
    for (index, path) in enumerate(paths)
        node = root
        for key in _canonical_update_path(model, _path_keys(path))
            node.index == 0 || throw(ArgumentError("Parameter paths overlap"))
            node = get!(node.children, key) do
                _ModelUpdateNode()
            end
        end
        node.index == 0 && isempty(node.children) || throw(ArgumentError("Parameter paths overlap or resolve to the same destination"))
        node.index = index
    end
    return root
end

function _apply_model_updates(object, node::_ModelUpdateNode, values; threaded_components::Bool=true, reference_model=nothing, tracking_steps=nothing)
    node.index != 0 && return values[node.index]
    changes = Dict{Symbol,Any}()
    changed = false
    entries = collect(node.children)
    replacements = Vector{Any}(undef, length(entries))
    function update_child(i)
        key, child = entries[i]
        replacements[i] = _apply_model_updates(_get(object, key), child, values; threaded_components, reference_model, tracking_steps)
    end
    if threaded_components && object isa AbstractVector{<:Component} && length(entries) > 1 && Threads.nthreads() > 1
        Threads.@threads for i in eachindex(entries)
            update_child(i)
        end
    else
        foreach(update_child, eachindex(entries))
    end
    for (i, (key, child)) in enumerate(entries)
        original = _get(object, key)
        replacement = replacements[i]
        changed |= replacement !== original
        # A selected computed property (gate duration) must retain its explicit
        # requested value when another input changes its inferred value.
        computed = !(object isa Union{NamedTuple,AbstractDict,AbstractArray,Tuple}) &&
            !hasfield(typeof(object), key)
        replacement !== original || computed || (object isa DeviceModel && key === :gates) || continue
        changes[key] = replacement
    end
    changed || return object
    return object isa DeviceModel ? _rebuild_model(object, changes; threaded_components, reference_model, tracking_steps) :
        _replace_properties(object, changes)
end

# Conservative fallback preserves specialized reconstruction hooks on custom
# objects. Components and gates specialize this to validate their final values.
function _replace_properties(object, changes)
    return foldl(pairs(changes); init=object) do candidate, (key, value)
        _replace_property(candidate, key, value)
    end
end

function _replace_properties(object::NamedTuple, changes)
    return merge(object, (; changes...))
end

function _replace_properties(object::AbstractDict, changes)
    result = copy(object)
    for (key, value) in changes
        result[_dict_key(object, key)] = value
    end
    return result
end

function _replace_properties(object::AbstractArray, changes)
    T = promote_type(eltype(object), (typeof(v) for v in values(changes))...)
    result = similar(object, T)
    copyto!(result, object)
    for (key, value) in changes
        result[_path_index(object, key)] = value
    end
    return result
end

_replace_properties(component::Component, changes) = _rebuild_component(component, changes)
