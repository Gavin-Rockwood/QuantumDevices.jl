"""
    ErfSquare(sigma;duration,ramp_time=2sigma)

Error-function edges and a constant plateau. Both physical ramps must fit duration.
Returns an Envelope with a physical reference duration and normalized shape payload.
"""
struct ErfSquare{T0,T1}
    sigma::T0
    ramp_time::T1
    ErfSquare(::Val{:raw},sigma,ramp_time)=new{typeof(sigma),typeof(ramp_time)}(sigma,ramp_time)
end
function ErfSquare(;duration=nothing,sigma,ramp_time=2sigma,rise=nothing,fall=nothing)
    _require_positive_width("duration",duration)
    _require_finite_real("sigma",sigma)
    _require_finite_real("ramp_time",ramp_time)
    sigma>0 || throw(ArgumentError("sigma must be positive"))
    0<=ramp_time<=duration/2 || throw(ArgumentError("Both ramps must fit duration"))
    return _wrap_shape(ErfSquare(Val(:raw),sigma/duration,ramp_time/duration),duration,rise,fall)
end
ErfSquare(x;duration=nothing,kwargs...)=ErfSquare(;duration,sigma=x,kwargs...)
_physical_parameters(s::ErfSquare,d)=(sigma=s.sigma*d,ramp_time=s.ramp_time*d,)
_shape_from_physical(::ErfSquare,p,d)=ErfSquare(;duration=d,sigma=p.sigma,ramp_time=p.ramp_time).shape
