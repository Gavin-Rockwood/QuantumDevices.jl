"""A shape piece with inferred/explicit duration and finite real or complex scale."""
struct Segment{E<:AbstractEnvelope,D<:Real,S}
    shape::E
    duration::D
    scale::S
    function Segment(shape::AbstractEnvelope;duration=nothing,scale=1)
        duration===nothing && (duration=_envelope_duration(shape))
        duration===nothing && throw(ArgumentError("Segment duration is required"))
        _require_positive_width("segment duration",duration)
        _require_finite("segment scale",scale)
        e=_stretch_envelope(shape,duration)
        return new{typeof(e),typeof(duration),typeof(scale)}(e,duration,scale)
    end
end
_segment(s::Segment)=s
function _segment(s::NamedTuple)
    haskey(s,:shape) && all(k->k in (:shape,:duration,:scale),keys(s)) || throw(ArgumentError("Piece requires shape and only optional duration/scale"))
    return Segment(s.shape;duration=get(s,:duration,nothing),scale=get(s,:scale,1))
end
"""Concatenated shape pieces; duration is the checked sum and cannot be stretched as a whole."""
struct EnvelopeSequence{S<:Tuple} <: AbstractEnvelope
    segments::S
    function EnvelopeSequence(pieces...)
        isempty(pieces) && throw(ArgumentError("EnvelopeSequence needs pieces"))
        segments=map(_segment,pieces)
        foldl(_checked_time_add,(s.duration for s in segments))
        return new{typeof(segments)}(segments)
    end
end
_envelope_duration(e::EnvelopeSequence)=foldl(_checked_time_add,(s.duration for s in e.segments))
_envelope_duration(::AbstractEnvelope)=nothing
function _stretch_envelope(e::EnvelopeSequence,d)
    d==_envelope_duration(e) || throw(ArgumentError("Edit segment durations; an EnvelopeSequence cannot be stretched"))
    return e
end
validate_envelope(e::EnvelopeSequence,d)=(_stretch_envelope(e,d);nothing)
function envelope_value(e::EnvelopeSequence,t,d)
    start=zero(t)
    for (i,s) in enumerate(e.segments)
        stop=_checked_time_add(start,s.duration)
        if t<stop || i==length(e.segments)
            return s.scale*s.shape(t==stop ? s.duration : t-start,s.duration)
        end
        start=stop
    end
end
function envelope_tstops(e::EnvelopeSequence,d)
    result=Real[]
    start=zero(d)
    for s in e.segments
        append!(result,[_checked_time_add(start,t) for t in _validated_envelope_tstops(s.shape,s.duration)])
        start=_checked_time_add(start,s.duration)
        push!(result,start)
    end
    return sort!(unique(filter(t->0<t<d,result)))
end
function Base.getproperty(e::EnvelopeSequence,k::Symbol)
    k===:duration && return _envelope_duration(e)
    k in (:rise,:fall) && return nothing
    k in (:rising,:falling) && throw(ArgumentError("EnvelopeSequence has no $k region"))
    return getfield(e,k)
end
Base.propertynames(e::EnvelopeSequence,private::Bool=false)=(:segments,:duration,:rise,:fall,:rising,:falling)
