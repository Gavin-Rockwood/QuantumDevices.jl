# Built-in shape formulas bind physical scalar accessors once.
function _shape_accessors(shape,d,base,selected,duration)
    defaults=_physical_parameters(shape,d)
    return map(keys(defaults)) do key
        value=getproperty(defaults,key)
        path=_path_join(base,key)
        # All dimensional shape fields scale with the authoritative duration.
        dimensional=key in (:sigma,:beta,:width,:center,:ramp_time,:rise_time,:fall_time)
        _scalar_accessor(value,path,selected;ratio=dimensional ? (iszero(d) ? getproperty(_physical_parameters(shape,one(d)),key) : value/d) : nothing,duration)
    end
end
function _compile_shape(s,d,base,selected,duration)
    !any(path -> startswith(path,base*"/"),selected) && return (p,t) -> _shape_value(s,t/duration(p))
    fields=_physical_parameters(s,d)
    accessors=_shape_accessors(s,d,base,selected,duration)
    a(k)=accessors[findfirst(==(k),keys(fields))]
    if s isa Constant
        return (p,t) -> one(t)
    elseif s isa Union{Gaussian,GaussianZero,DRAG}
        sigma=a(:sigma)
        s isa Gaussian && return (p,t) -> exp(-(t-duration(p)/2)^2/(2sigma(p)^2))
        s isa GaussianZero && return (p,t) -> _zero_gaussian(sigma(p),t,duration(p))
        beta=a(:beta)
        return (p,t) -> complex(_zero_gaussian(sigma(p),t,duration(p)),beta(p)*_zero_gaussian_derivative(sigma(p),t,duration(p)))
    elseif s isa Union{SineSquared,Cosine,Blackman}
        return (p,t) -> _shape_value(s,t/duration(p))
    elseif s isa Union{Sech,Bump}
        width,center=a(:width),a(:center)
        s isa Sech && return (p,t) -> inv(cosh((t-center(p))/width(p)))
        k=a(:k)
        return (p,t) -> begin
            x=(t-center(p))/width(p)
            abs(x)>=1 ? zero(float(x)) : exp(k(p)*x^2/(x^2-1))
        end
    elseif s isa Union{GaussianSquare,ErfSquare}
        sigma,ramp=a(:sigma),a(:ramp_time)
        return (p,t) -> begin
            r=ramp(p); distance=min(t,duration(p)-t)
            (iszero(r) || distance>=r) && return one(float(t))
            s isa GaussianSquare && return _zero_gaussian(sigma(p),distance,2r)
            distance<=0 && return zero(float(t))
            h=r/(2sqrt(2)*sigma(p)); denominator=erf(h)
            iszero(denominator) ? distance/r : (erf((distance-r/2)/(sqrt(2)*sigma(p)))+denominator)/(2denominator)
        end
    elseif s isa Region
        lo,hi=s.interval
        source_duration=p -> duration(p)/(hi-lo)
        source=_compile_shape(s.shape,d/(hi-lo),base,selected,source_duration)
        return (p,t) -> source(p,t+lo*source_duration(p))
    elseif s isa Slepian
        any(path -> startswith(path,base*"/time_bandwidth"),selected) && throw(ArgumentError("$base/time_bandwidth requires rebuilding Slepian weights"))
        return (p,t) -> _shape_value(s,t/duration(p))
    elseif s isa InterpolatedEnvelope
        times=map(Tuple(eachindex(s.times))) do i
            _scalar_accessor(s.times[i]*d,_path_join(base,"times/$i"),selected;ratio=s.times[i],duration)
        end
        values=map(Tuple(eachindex(s.values))) do i
            _scalar_accessor(s.values[i],_path_join(base,"values/$i"),selected)
        end
        return (p,t) -> begin
            t<=0 && return values[1](p)
            t>=duration(p) && return values[end](p)
            low,high=1,length(times)
            while low+1<high
                middle=(low+high)÷2
                if times[middle](p)<=t; low=middle; else; high=middle; end
            end
            x,y=times[low](p),times[low+1](p)
            v,w=values[low](p),values[low+1](p)
            v+(w-v)*(t-x)/(y-x)
        end
    elseif s isa RampedFlattop
        rise,fall=a(:rise_time),a(:fall_time)
        up=_compile_ramp(s.ramp_up,s.rise_time*d,_path_join(base,:ramp_up),selected,rise,true)
        down=_compile_ramp(s.ramp_down,s.fall_time*d,_path_join(base,:ramp_down),selected,fall,false)
        return (p,t) -> begin
            r,f=rise(p),fall(p); total=duration(p)
            t<r && return up(p,t)
            (t<total-f || iszero(f)) && return one(t)
            down(p,t-(total-f))
        end
    elseif s isa CustomShape
        # Explicit custom state is immutable; its callback owns evaluation semantics.
        any(path -> startswith(path,base*"/parameters/"),selected) && throw(ArgumentError("runtime custom-state fields are not yet supported: $base"))
        return (p,t) -> s.callable(s.parameters,t/duration(p))
    end
    throw(ArgumentError("no scalar waveform compiler for $(typeof(s))"))
end
function _compile_ramp(s,d,base,selected,duration,rising)
    s isa NamedTuple && return _compile_shape(s.shape,d,base,selected,duration)
    return (p,t) -> begin
        u=t/duration(p); x=rising ? u : 1-u
        s===:linear ? x : s===:smoothstep ? x*x*(3-2x) : sinpi(x/2)^2
    end
end
function _compile_envelope(e::Envelope,base,selected,duration)
    return _compile_shape(e.shape,e.duration,base,selected,duration)
end
function _compile_envelope(e::EnvelopeSequence,base,selected,duration)
    segments=map(Tuple(enumerate(e.segments))) do (i,s)
        prefix=_path_join(base,"segments/$i")
        d=s.shape isa EnvelopeSequence ? _compiled_envelope_duration(s.shape,_path_join(prefix,:shape),selected) : _scalar_accessor(s.duration,_path_join(prefix,:duration),selected)
        f=_compile_envelope(s.shape,_path_join(prefix,:shape),selected,d)
        scale=_scalar_accessor(s.scale,_path_join(prefix,:scale),selected)
        (d,f,scale)
    end
    return (p,t) -> begin
        start=zero(t)
        for (i,(d,f,scale)) in enumerate(segments)
            stop=start+d(p)
            (t<stop || i==length(segments)) && return scale(p)*f(p,t-start)
            start=stop
        end
    end
end
function _compile_envelope(e::LegacyEnvelope,base,selected,duration)
    any(path -> startswith(path,base*"/"),selected) && throw(ArgumentError("legacy physical-envelope fields require candidate rebuilding: $base"))
    source=_legacy_value(e)
    return (p,t) -> _LegacyPhysical.envelope_value(source,t,duration(p))
end


# Extensions keep their own pure envelope_value implementation for fixed state.
function _compile_envelope(e::AbstractEnvelope,base,selected,duration)
    any(path -> startswith(path,base*"/"),selected) && throw(ArgumentError("custom envelope fields require rebuilding: $base"))
    return (p,t) -> envelope_value(e,t,duration(p))
end
