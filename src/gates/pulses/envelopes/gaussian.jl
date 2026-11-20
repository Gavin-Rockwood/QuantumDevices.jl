"""
    Gaussian(;duration,sigma)

Gaussian centered at half its reference duration. sigma is physical time; stretch scales it.
Returns an Envelope with a physical reference duration and normalized shape payload.
"""
struct Gaussian{T0}
    sigma::T0
    Gaussian(::Val{:raw},sigma)=new{typeof(sigma)}(sigma)
end
function Gaussian(;duration=nothing,sigma,rise=_default_rise(duration),fall=_default_fall(duration))
    _require_positive_width("duration",duration)
    _require_finite_real("sigma",sigma)
    sigma>0 || throw(ArgumentError("sigma must be positive"))
    return _wrap_shape(Gaussian(Val(:raw),sigma/duration),duration,rise,fall)
end
function Gaussian(x;duration=nothing,kwargs...)
    haskey(kwargs,:center) && throw(ArgumentError("Gaussian is centered at half its duration; use Pulse.delay for placement or an explicit custom shape for an off-center legacy Gaussian"))
    return Gaussian(;duration,sigma=x,kwargs...)
end
_physical_parameters(s::Gaussian,d)=(sigma=s.sigma*d,)
_shape_from_physical(::Gaussian,p,d)=Gaussian(;duration=d,sigma=p.sigma).shape
