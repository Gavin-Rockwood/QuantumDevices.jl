"""
    Sech(width;duration,center=duration/2)

Hyperbolic secant envelope with physical width and center.
Returns an Envelope with a physical reference duration and normalized shape payload.
"""
struct Sech{T0,T1}
    width::T0
    center::T1
    Sech(::Val{:raw},width,center)=new{typeof(width),typeof(center)}(width,center)
end
function Sech(;duration=nothing,width,center=nothing,rise=nothing,fall=nothing)
    _require_positive_width("duration",duration)
    center===nothing && (center=duration/2)
    _require_finite_real("width",width)
    _require_finite_real("center",center)
    width>0 || throw(ArgumentError("width must be positive"))
    0<=center<=duration || throw(ArgumentError("center must lie inside duration"))
    return _wrap_shape(Sech(Val(:raw),width/duration,center/duration),duration,rise,fall)
end
Sech(x;duration=nothing,kwargs...)=Sech(;duration,width=x,kwargs...)
_physical_parameters(s::Sech,d)=(width=s.width*d,center=s.center*d,)
_shape_from_physical(::Sech,p,d)=Sech(;duration=d,width=p.width,center=p.center).shape
