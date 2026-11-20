"""
    GaussianSquare(sigma;duration,ramp_time=2sigma)

Gaussian edges and a constant plateau. Both physical ramps must fit duration.
Returns an Envelope with a physical reference duration and normalized shape payload.
"""
struct GaussianSquare{T0,T1}
    sigma::T0
    ramp_time::T1
    GaussianSquare(::Val{:raw},sigma,ramp_time)=new{typeof(sigma),typeof(ramp_time)}(sigma,ramp_time)
end
function GaussianSquare(;duration=nothing,sigma,ramp_time=2sigma,rise=nothing,fall=nothing)
    _require_positive_width("duration",duration)
    _require_finite_real("sigma",sigma)
    _require_finite_real("ramp_time",ramp_time)
    sigma>0 || throw(ArgumentError("sigma must be positive"))
    0<=ramp_time<=duration/2 || throw(ArgumentError("Both ramps must fit duration"))
    return _wrap_shape(GaussianSquare(Val(:raw),sigma/duration,ramp_time/duration),duration,rise,fall)
end
GaussianSquare(x;duration=nothing,kwargs...)=GaussianSquare(;duration,sigma=x,kwargs...)
_physical_parameters(s::GaussianSquare,d)=(sigma=s.sigma*d,ramp_time=s.ramp_time*d,)
_shape_from_physical(::GaussianSquare,p,d)=GaussianSquare(;duration=d,sigma=p.sigma,ramp_time=p.ramp_time).shape
