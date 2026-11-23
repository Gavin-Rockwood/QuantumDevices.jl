const _RAMP_NAMES=(:sine_squared,:linear,:smoothstep)
# Source envelopes contribute normalized payloads, not a second duration owner.
_ramp_source(s::NamedTuple)=s
_ramp_source(s::Envelope)=(shape=s.shape,rise=s.rise,fall=s.fall)
function _ramp_source(s::Symbol)
    s in _RAMP_NAMES || throw(ArgumentError("Unknown ramp"))
    return s
end
function _ramp_value(s,u,rising)
    s isa NamedTuple && return _shape_value(s.shape,u)
    x=rising ? u : 1-u
    s===:linear && return x
    s===:smoothstep && return x*x*(3-2x)
    return sinpi(x/2)^2
end
"""
    RampedFlattop(ramp_time;duration,ramp=:sine_squared)

Direct normalized rise/plateau/fall envelope. Envelope ramp sources sample forward; supply rising/falling views explicitly.
Returns an Envelope with a physical reference duration and normalized shape payload.
"""
struct RampedFlattop{R,F,U,D}
    rise_time::R
    fall_time::F
    ramp_up::U
    ramp_down::D
    RampedFlattop(::Val{:raw},r,f,u,d)=new{typeof(r),typeof(f),typeof(u),typeof(d)}(r,f,u,d)
end
function RampedFlattop(r;duration=nothing,ramp=:sine_squared,rise_time=r,fall_time=r,ramp_up=ramp,ramp_down=ramp,rise=nothing,fall=nothing)
    _require_positive_width("duration",duration)
    _require_finite_real("rise_time",rise_time)>=0 || throw(ArgumentError("Negative rise"))
    _require_finite_real("fall_time",fall_time)>=0 || throw(ArgumentError("Negative fall"))
    _checked_time_add(rise_time,fall_time)<=duration || throw(ArgumentError("Ramps must fit duration"))
    s=RampedFlattop(Val(:raw),rise_time/duration,fall_time/duration,_ramp_source(ramp_up),_ramp_source(ramp_down))
    return _wrap_shape(s,duration,rise,fall)
end
function _shape_value(s::RampedFlattop,u)
    u<s.rise_time && return _ramp_value(s.ramp_up,u/s.rise_time,true)
    (u<1-s.fall_time || iszero(s.fall_time)) && return one(u)
    return _ramp_value(s.ramp_down,u==1 ? one(u) : (u-(1-s.fall_time))/s.fall_time,false)
end
function _shape_tstops(s::RampedFlattop)
    result=Real[]
    s.rise_time>0 && push!(result,s.rise_time)
    s.fall_time>0 && push!(result,1-s.fall_time)
    s.ramp_up isa NamedTuple && append!(result,s.rise_time.*collect(_shape_tstops(s.ramp_up.shape)))
    s.ramp_down isa NamedTuple && append!(result,1-s.fall_time .+ s.fall_time.*collect(_shape_tstops(s.ramp_down.shape)))
    return sort!(unique(result))
end
_ramp_physical(s::Symbol,d)=s
_ramp_physical(s::NamedTuple,d)=Envelope(Val(:raw),s.shape,d,s.rise,s.fall)
_physical_parameters(s::RampedFlattop,d)=(rise_time=s.rise_time*d,fall_time=s.fall_time*d,
    ramp_up=_ramp_physical(s.ramp_up,s.rise_time*d),
    ramp_down=_ramp_physical(s.ramp_down,s.fall_time*d))
_shape_from_physical(::RampedFlattop,p,d)=RampedFlattop(p.rise_time;duration=d,fall_time=p.fall_time,ramp_up=p.ramp_up,ramp_down=p.ramp_down).shape
