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
Built-in components rebuild their basis, and model physical-input updates rebuild
operators, the Hamiltonian, and dressed states while retaining named gates.
Derived model/component fields cannot be directly replaced.

```jldoctest
julia> setpath((amplitude=0.2,), "amplitude", 0.4)
(amplitude = 0.4,)
```
"""
setpath(object, path::AbstractString, value) = _setpath(object, _path_keys(path), value)

"""
    setpath!(object, path::AbstractString, value)

Update a mutable root (dictionary, array, or mutable struct) and return it.
Immutable descendants are reconstructed first. An immutable root such as a
[`DeviceGate`](@ref) or [`DeviceModel`](@ref) requires [`setpath`](@ref).
Path semantics and reconstruction restrictions are the same as `setpath`.
"""
function setpath!(object, path::AbstractString, value)
    (object isa Union{AbstractArray,AbstractDict} || ismutabletype(typeof(object))) ||
        throw(ArgumentError("Cannot mutate an immutable root; use object = setpath(object, path, value)"))
    keys = _path_keys(path)
    key = first(keys)
    if !(object isa Union{AbstractArray,AbstractDict}) && !hasproperty(object, key)
        hasproperty(object, :parameters) || throw(KeyError(key))
        replacement = _setpath(object.parameters, keys, value)
        key = :parameters
    else
        replacement = length(keys) == 1 ? value : _setpath(_get(object, key), keys[2:end], value)
    end
    if object isa AbstractDict
        object[_dict_key(object, key; insert = true)] = replacement
    elseif object isa AbstractArray
        object[_path_index(object, key)] = replacement
    else
        candidate = _replace_property(object, key, replacement)
        # Convert before mutation; failed conversions leave the root unchanged.
        converted = convert(fieldtype(typeof(object), key), getfield(candidate, key))
        setproperty!(object, key, converted)
    end
    return object
end
