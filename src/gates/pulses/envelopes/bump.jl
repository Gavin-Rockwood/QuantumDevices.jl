"""
    Bump(width;duration,k=2,center=duration/2)

Compact support envelope with physical half-width and dimensionless sharpness k.
Returns an Envelope with a physical reference duration and normalized shape payload.
"""
struct Bump{T0,T1,T2}
    width::T0
    k::T1
    center::T2
    Bump(::Val{:raw},width,k,center)=new{typeof(width),typeof(k),typeof(center)}(width,k,center)
end
function Bump(;duration=nothing,width,k=2,center=nothing,rise=nothing,fall=nothing)
    _require_positive_width("duration",duration)
    center===nothing && (center=duration/2)
    _require_finite_real("width",width)
    _require_finite_real("k",k)
    _require_finite_real("center",center)
    width>0 || throw(ArgumentError("width must be positive"))
    k>0 || throw(ArgumentError("k must be positive"))
    0<=center<=duration || throw(ArgumentError("center must lie inside duration"))
    return _wrap_shape(Bump(Val(:raw),width/duration,k,center/duration),duration,rise,fall)
end
Bump(x;duration=nothing,kwargs...)=Bump(;duration,width=x,kwargs...)
_physical_parameters(s::Bump,d)=(width=s.width*d,k=s.k,center=s.center*d,)
_shape_from_physical(::Bump,p,d)=Bump(;duration=d,width=p.width,k=p.k,center=p.center).shape
