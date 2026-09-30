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
