"""
    DRAG(sigma;duration,beta=0)

Gaussian-zero I with derivative Q. Physical beta and sigma stretch with duration.
Returns an Envelope with a physical reference duration and normalized shape payload.
"""
struct DRAG{T0,T1}
    sigma::T0
    beta::T1
    DRAG(::Val{:raw},sigma,beta)=new{typeof(sigma),typeof(beta)}(sigma,beta)
end
function DRAG(;duration=nothing,sigma,beta=0,rise=nothing,fall=nothing)
    _require_positive_width("duration",duration)
    _require_finite_real("sigma",sigma)
    _require_finite_real("beta",beta)
    sigma>0 || throw(ArgumentError("sigma must be positive"))
    return _wrap_shape(DRAG(Val(:raw),sigma/duration,beta/duration),duration,rise,fall)
end
DRAG(x;duration=nothing,kwargs...)=DRAG(;duration,sigma=x,kwargs...)
_physical_parameters(s::DRAG,d)=(sigma=s.sigma*d,beta=s.beta*d,)
_shape_from_physical(::DRAG,p,d)=DRAG(;duration=d,sigma=p.sigma,beta=p.beta).shape
