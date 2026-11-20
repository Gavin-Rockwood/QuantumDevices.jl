"""
    GaussianZero(sigma;duration)

Centered Gaussian with zero endpoints and unit peak.
Returns an Envelope with a physical reference duration and normalized shape payload.
"""
struct GaussianZero{T0}
    sigma::T0
    GaussianZero(::Val{:raw},sigma)=new{typeof(sigma)}(sigma)
end
function GaussianZero(;duration=nothing,sigma,rise=nothing,fall=nothing)
    _require_positive_width("duration",duration)
    _require_finite_real("sigma",sigma)
    sigma>0 || throw(ArgumentError("sigma must be positive"))
    return _wrap_shape(GaussianZero(Val(:raw),sigma/duration),duration,rise,fall)
end
GaussianZero(x;duration=nothing,kwargs...)=GaussianZero(;duration,sigma=x,kwargs...)
_physical_parameters(s::GaussianZero,d)=(sigma=s.sigma*d,)
_shape_from_physical(::GaussianZero,p,d)=GaussianZero(;duration=d,sigma=p.sigma).shape
