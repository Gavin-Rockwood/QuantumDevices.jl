"""A captured static conversion from ordered lab voltages to absolute normalized flux."""
abstract type AbstractControlMap end

# Control maps own immutable numeric tables and supported frozen metadata.
_map_real(x) = x isa Real && !(x isa Union{Bool,Num}) && isfinite(x)
function _map_numbers(values, label)
    all(_map_real, values) || throw(ArgumentError("$label must contain finite real numbers"))
    return promote(values...)
end
function _map_axes(inputs, outputs)
    ins, outs = Tuple(inputs), Tuple(outputs)
    !isempty(ins) && all(x -> x isa Symbol, ins) && allunique(ins) || throw(ArgumentError("map inputs must be distinct Symbols"))
    !isempty(outs) && all(is_parameter, outs) && allunique(outs) || throw(ArgumentError("map outputs must be distinct atomic symbolic variables"))
    return ins, outs
end
function _map_bounds(bounds, count, label)
    isnothing(bounds) && return nothing
    length(bounds) == count || throw(ArgumentError("$label width differs from map axis"))
    return Tuple(map(bounds) do pair
        length(pair) == 2 || throw(ArgumentError("$label requires (lower,upper) pairs"))
        lo, hi = _map_numbers(Tuple(pair), label)
        lo <= hi || throw(ArgumentError("$label lower bound exceeds upper bound"))
        (lo, hi)
    end)
end
include("AffineControlMap.jl")
include("TabulatedControlMap.jl")
include("evaluation.jl")

_control_map_content(map::AbstractControlMap) = Tuple(getfield(map, k) for k in fieldnames(typeof(map)))
_control_map_derivation_content(value) = value
_control_map_content(map::TabulatedControlMap) = Tuple(k===:derivation ? _control_map_derivation_content(map.derivation) : getfield(map,k) for k in fieldnames(typeof(map)))
Base.isequal(a::AbstractControlMap, b::AbstractControlMap) = typeof(a).name === typeof(b).name && isequal(_control_map_content(a), _control_map_content(b))
Base.:(==)(a::AbstractControlMap, b::AbstractControlMap) = isequal(a,b)
Base.hash(map::AbstractControlMap,h::UInt) = hash(_control_map_content(map),hash(nameof(typeof(map)),h))
_copy_control_map(map::Union{AffineControlMap,TabulatedControlMap}) = map
_copy_control_map(map::AbstractControlMap) = throw(ArgumentError("$(typeof(map)) requires an explicit owned-copy implementation"))

# Membership is symbolic equality, independent of component naming conventions.
function _validate_control_map_outputs(map, variables; context="control map $(map.name)")
    all(output -> any(isequal(output), variables), map.outputs) || throw(ArgumentError("$context has outputs outside selected physical variables"))
    return nothing
end
