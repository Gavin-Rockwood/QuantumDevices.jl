# Only operator declaration maps have this lifecycle: mutable draft, frozen registration.
# Tuple storage drops all dictionary aliases when the declaration is registered.
mutable struct OperatorDeclarations{V} <: AbstractDict{String,V}
    entries::Union{Dict{String,V},Tuple}
end
OperatorDeclarations(source::AbstractDict{String,V}) where V = OperatorDeclarations{V}(Dict(source))
Base.length(map::OperatorDeclarations) = length(map.entries)
Base.iterate(map::OperatorDeclarations, state...) = iterate(map.entries, state...)
function Base.getindex(map::OperatorDeclarations, key)
    for (name,value) in map
        isequal(name,key) && return value
    end
    throw(KeyError(key))
end
Base.haskey(map::OperatorDeclarations, key) = any(p -> isequal(first(p),key), map)
function editable_operator_entries(map::OperatorDeclarations)
    map.entries isa Dict || throw(ArgumentError("registered operator declarations are immutable"))
    return map.entries
end
Base.setindex!(map::OperatorDeclarations, value, key) = setindex!(editable_operator_entries(map),value,key)
# A factory may be replaced for an existing declaration; its local key is fixed.
# Physical identity is carried by the separate, fully frozen symbolic map.
function Base.setindex!(map::OperatorDeclarations{Function}, value::Function, key::String)
    applicable(value, 1) || throw(ArgumentError("operator generators must accept a dimension"))
    if map.entries isa Tuple
        haskey(map,key) || throw(ArgumentError("registered operator generator keys are immutable"))
        setfield!(map,:entries,Tuple(name => (name == key ? value : generator) for (name,generator) in map))
        return value
    end
    return setindex!(map.entries,value,key)
end
Base.delete!(map::OperatorDeclarations, key) = (delete!(editable_operator_entries(map),key); map)
Base.empty!(map::OperatorDeclarations) = (empty!(editable_operator_entries(map)); map)
freeze_operator_declarations!(map::OperatorDeclarations) = (map.entries = Tuple(pairs(map)); map)
function Base.setproperty!(map::OperatorDeclarations, name::Symbol, value)
    map.entries isa Tuple && throw(ArgumentError("registered operator declarations are immutable"))
    return setfield!(map,name,value)
end
