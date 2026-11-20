# Regions share the source payload; nested selections flatten to one interval.
struct Region{S,I}
    shape::S
    interval::I
end
function _region(e::Envelope,interval)
    a,b=interval
    source=e.shape
    if source isa Region
        x,y=source.interval
        a,b=x+(y-x)*a,x+(y-x)*b
        source=source.shape
    end
    x,y=interval
    remap(r)=r===nothing || max(r[1],x)>=min(r[2],y) ? nothing : ((max(r[1],x)-x)/(y-x),(min(r[2],y)-x)/(y-x))
    return _raw_envelope(Region(source,(a,b)),_checked_time_mul(y-x,e.duration),remap(e.rise),remap(e.fall))
end
_shape_value(s::Region,u)=_shape_value(s.shape,s.interval[1]+(s.interval[2]-s.interval[1])*u)
_shape_tstops(s::Region)=Tuple((u-s.interval[1])/(s.interval[2]-s.interval[1]) for u in _shape_tstops(s.shape) if s.interval[1]<u<s.interval[2])
_physical_parameters(s::Region,d)=_physical_parameters(s.shape,d/(s.interval[2]-s.interval[1]))
_shape_from_physical(s::Region,p,d)=Region(_shape_from_physical(s.shape,p,d/(s.interval[2]-s.interval[1])),s.interval)
function Base.getproperty(e::Envelope,key::Symbol)
    key in (:shape,:duration,:rise,:fall) && return getfield(e,key)
    if key in (:rising,:falling)
        interval=getfield(e,key===:rising ? :rise : :fall)
        interval===nothing && throw(ArgumentError("Envelope has no $key region"))
        return _region(e,interval)
    end
    p=_physical_parameters(getfield(e,:shape),getfield(e,:duration))
    haskey(p,key) || throw(ErrorException("Envelope has no property $key"))
    return getproperty(p,key)
end
Base.propertynames(e::Envelope,private::Bool=false)=(:shape,:duration,:rise,:fall,:rising,:falling,keys(_physical_parameters(e.shape,e.duration))...)
function _validate_shape(s::Region)
    _normalized_interval(s.interval)
    _validate_shape(s.shape)
    return nothing
end
