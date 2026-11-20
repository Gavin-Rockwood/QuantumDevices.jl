_default_rise(d)=d===nothing ? nothing : (zero(d),d/2)
_default_fall(d)=d===nothing ? nothing : (d/2,d)
"""Normalized shape with a physical reference duration and optional region metadata."""
struct Envelope{S,D,R,F} <: AbstractEnvelope
    shape::S
    duration::D
    rise::R
    fall::F
    Envelope(::Val{:raw}, s::S,d::D,r::R,f::F) where {S,D,R,F} = new{S,D,R,F}(s,d,r,f)
end
struct CustomShape{F,P<:NamedTuple}
    callable::F
    parameters::P
end
function _normalized_interval(interval)
    interval === nothing && return nothing
    length(interval)==2 || throw(ArgumentError("Region needs two endpoints"))
    a,b=interval
    all(x->x isa Real && isfinite(x),(a,b)) && 0<=a<b<=1 || throw(ArgumentError("Invalid normalized region"))
    return (a,b)
end
function _raw_envelope(shape,duration,rise=nothing,fall=nothing)
    _require_positive_width("duration",duration)
    e=Envelope(Val(:raw),shape,duration,_normalized_interval(rise),_normalized_interval(fall))
    validate_envelope(e,duration)
    return e
end
function _physical_interval(interval,duration)
    interval===nothing && return nothing
    length(interval)==2 || throw(ArgumentError("Region needs two endpoints"))
    return _normalized_interval(Tuple(x/duration for x in interval))
end
function _wrap_shape(shape,duration,rise=nothing,fall=nothing)
    _require_positive_width("duration",duration)
    return _raw_envelope(shape,duration,_physical_interval(rise,duration),_physical_interval(fall,duration))
end
function Envelope(f,p::NamedTuple=(;);duration=nothing,rise=nothing,fall=nothing)
    applicable(f,p,0.0) || throw(ArgumentError("Custom Envelope requires f(parameters,u), with normalized u; explicitly adapt old f(p,t,d) callbacks"))
    return _wrap_shape(CustomShape(f,p),duration,rise,fall)
end
function _validate_shape(s::CustomShape)
    applicable(s.callable,s.parameters,0.0) || throw(ArgumentError("Custom shape requires f(parameters,u)"))
    return nothing
end
_shape_value(s::CustomShape,u)=s.callable(s.parameters,u)
_physical_parameters(s::CustomShape,d)=(callable=s.callable,parameters=s.parameters)
_shape_from_physical(::CustomShape,p,d)=CustomShape(p.callable,p.parameters)
_shape_tstops(s)=()
_physical_parameters(s,d)=(;)
_shape_from_physical(s,p,d)=s
function validate_envelope(e::Envelope,d)
    _require_positive_width("duration",d)
    _validate_shape(e.shape)
    _require_finite("Envelope value",_shape_value(e.shape,0.0))
    _require_finite("Envelope value",_shape_value(e.shape,1.0))
    return nothing
end
function envelope_value(e::Envelope,t,d)
    _require_positive_width("evaluation duration",d)
    return _shape_value(e.shape,t/d)
end
envelope_tstops(e::Envelope,d)=[_checked_time_mul(u,d) for u in _shape_tstops(e.shape)]
_envelope_duration(e::Envelope)=getfield(e,:duration)
_stretch_envelope(e::Envelope,d)=_raw_envelope(e.shape,d,e.rise,e.fall)
