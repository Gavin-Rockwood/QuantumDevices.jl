"""
    InterpolatedEnvelope(times,values;duration,interpolation=:linear)

Linear real/complex samples from physical time zero to duration. Input nodes normalize
once; duration-only changes stretch nodes. Physical times/values paths rebuild atomically.
"""
struct InterpolatedEnvelope{T,V}
    times::T
    values::V
    interpolation::Symbol
    InterpolatedEnvelope(::Val{:raw},t,v,i)=new{typeof(t),typeof(v)}(t,v,i)
end
function InterpolatedEnvelope(times,values;duration=nothing,interpolation=:linear,rise=nothing,fall=nothing)
    _require_positive_width("duration",duration)
    length(times)==length(values)>=2 || throw(ArgumentError("Need matching samples"))
    all(t->t isa Real && isfinite(t),times) && all(isfinite,values) || throw(ArgumentError("Nonfinite samples"))
    first(times)==0 && last(times)==duration && all(>(0),diff(collect(times))) || throw(ArgumentError("Samples must increase from zero to duration"))
    interpolation===:linear || throw(ArgumentError("Only linear interpolation is supported"))
    return _wrap_shape(InterpolatedEnvelope(Val(:raw),Tuple(t/duration for t in times),Tuple(values),interpolation),duration,rise,fall)
end
function _shape_value(s::InterpolatedEnvelope,u)
    u<=0 && return first(s.values)
    u>=1 && return last(s.values)
    low,high=1,length(s.times)
    while low+1<high
        middle=(low+high)÷2
        if s.times[middle]<=u
            low=middle
        else
            high=middle
        end
    end
    i=low
    return s.values[i]+(s.values[i+1]-s.values[i])*(u-s.times[i])/(s.times[i+1]-s.times[i])
end
_shape_tstops(s::InterpolatedEnvelope)=s.times[2:end-1]
_physical_parameters(s::InterpolatedEnvelope,d)=(times=collect(s.times).*d,values=collect(s.values),interpolation=s.interpolation)
_shape_from_physical(::InterpolatedEnvelope,p,d)=InterpolatedEnvelope(p.times,p.values;duration=d,interpolation=p.interpolation).shape
