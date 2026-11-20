# Recipe-owned collections use immutable tuples. This storage is confined to
# snapshots and provenance; editable drafts retain ordinary dictionaries/vectors.
struct RecipeVector{T, S <: Tuple} <: AbstractVector{T}
    items::S
end
# Freeze vector elements in their original order.
function RecipeVector(items)
    values = Tuple(freeze_recipe_value(item) for item in items)
    T = isempty(values) ? eltype(items) : reduce(typejoin, typeof.(values))
    return RecipeVector{T, typeof(values)}(values)
end

struct RecipeArray{T, N, S <: Tuple} <: AbstractArray{T, N}
    items::S
    dimensions::NTuple{N, Int}
end
# Preserve array shape while freezing its linear element storage.
function RecipeArray(source::AbstractArray{T, N}) where {T, N}
    items = Tuple(freeze_recipe_value(item) for item in source)
    Value = isempty(items) ? T : reduce(typejoin, typeof.(items))
    return RecipeArray{Value, N, typeof(items)}(items, size(source))
end

struct RecipeMap{K, V, S <: Tuple} <: AbstractDict{K, V}
    items::S
end
# Preserve dictionary iteration order and freeze keys and values.
function RecipeMap(source::AbstractDict{K, V}) where {K, V}
    items = Tuple(freeze_recipe_value(key) => freeze_recipe_value(value) for (key, value) in source)
    Value = isempty(items) ? V : reduce(typejoin, typeof(last(item)) for item in items)
    return RecipeMap{K, Value, typeof(items)}(items)
end


# Supported collections are copied recursively; callables retain their identity.
freeze_recipe_value(value::AbstractDict) = RecipeMap(value)
freeze_recipe_value(value::AbstractVector) = RecipeVector(value)
freeze_recipe_value(value::AbstractArray) = RecipeArray(value)
freeze_recipe_value(value::Tuple) = map(freeze_recipe_value, value)
freeze_recipe_value(value::Pair) =
    freeze_recipe_value(first(value)) => freeze_recipe_value(last(value))
freeze_recipe_value(value::NamedTuple) = map(freeze_recipe_value, value)
freeze_recipe_value(value::Function) = value
recipe_value_stable(value::Union{Function, Number, Symbol, String, OperatorGenerator}) = true
# Immutable custom values are accepted only when their fields are also stable.
function recipe_value_stable(value)
    ismutable(value) && return false
    return all(recipe_value_stable(getfield(value, field)) for field = 1:fieldcount(typeof(value)))
end
# Reject unsupported mutable payloads rather than silently sharing them.
function freeze_recipe_value(value)
    recipe_value_stable(value) || throw(
        ArgumentError(
            "unsupported mutable recipe value $(typeof(value)); use an immutable value or supported collection",
        ),
    )
    return value
end
# Drafts recover ordinary mutable collections from captured storage.
thaw_recipe_value(value::RecipeMap) = Dict(key => thaw_recipe_value(item) for (key, item) in value)
thaw_recipe_value(value::RecipeVector) = [thaw_recipe_value(item) for item in value]
thaw_recipe_value(value::Tuple) = map(thaw_recipe_value, value)
thaw_recipe_value(value::NamedTuple) = map(thaw_recipe_value, value)
thaw_recipe_value(value) = value

thaw_recipe_value(value::Pair) = thaw_recipe_value(first(value)) => thaw_recipe_value(last(value))

# Recover multidimensional array shape from linear snapshot storage.
function thaw_recipe_value(value::RecipeArray)
    return reshape([thaw_recipe_value(item) for item in value.items], size(value))
end

# Immutable recipe snapshot collection access
Base.size(items::RecipeVector) = (length(items.items),)
Base.IndexStyle(::Type{<:RecipeVector}) = IndexLinear()
Base.getindex(items::RecipeVector, index::Int) = items.items[index]
Base.copy(items::RecipeVector{T}) where {T} = T[items...]
Base.similar(::RecipeVector, ::Type{T}, dimensions::Dims) where {T} = Array{T}(undef, dimensions)
Base.length(map::RecipeMap) = length(map.items)
Base.iterate(map::RecipeMap, state...) = iterate(map.items, state...)
Base.haskey(map::RecipeMap, key) = any(isequal(first(item), key) for item in map.items)
Base.get(map::RecipeMap, key, default) = haskey(map, key) ? map[key] : default
Base.copy(map::RecipeMap) = Dict(map)
Base.merge(map::RecipeMap, other::AbstractDict) = merge(Dict(map), other)
function Base.getindex(map::RecipeMap, key)
    for (stored, value) in map
        isequal(stored, key) && return value
    end
    return throw(KeyError(key))
end

Base.show(io::IO, map::RecipeMap) = show(io, Dict(map))
Base.show(io::IO, items::RecipeVector) = show(io, collect(items))

Base.size(items::RecipeArray) = items.dimensions
Base.IndexStyle(::Type{<:RecipeArray}) = IndexLinear()
Base.getindex(items::RecipeArray, index::Int) = items.items[index]
Base.similar(::RecipeArray, ::Type{T}, dimensions::Dims) where {T} = Array{T}(undef, dimensions)
