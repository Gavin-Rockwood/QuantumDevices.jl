# Scalar accessors are assembled once; callbacks only read literal runtime keys.
_constant_accessor(x) = p -> x
function _scalar_accessor(value,path,selected;ratio=nothing,duration=nothing)
    path in selected && return p -> p[path]
    ratio === nothing && return _constant_accessor(value)
    return p -> ratio*duration(p)
end
_path_join(a,b) = isempty(a) ? string(b) : "$a/$b"
function _compiled_duration(r::PulseRecipe,base,selected)
    r.envelope isa EnvelopeSequence && return _compiled_envelope_duration(r.envelope,_path_join(base,:envelope),selected)
    return _scalar_accessor(r.duration,_path_join(base,:duration),selected)
end
function _compiled_envelope_duration(e,base,selected)
    ds=map(Tuple(enumerate(e.segments))) do (i,s)
        prefix=_path_join(base,"segments/$i")
        s.shape isa EnvelopeSequence ? _compiled_envelope_duration(s.shape,_path_join(prefix,:shape),selected) : _scalar_accessor(s.duration,_path_join(prefix,:duration),selected)
    end
    return p -> sum(f(p) for f in ds)
end
function _compiled_duration(r::PulseSequenceRecipe,base,selected)
    ds=map(Tuple(enumerate(r.pulses))) do (i,child)
        prefix=_path_join(base,"pulses/$i")
        d=_compiled_duration(child,prefix,selected)
        delay=_scalar_accessor(child.delay,_path_join(prefix,:delay),selected)
        p -> d(p)+delay(p)
    end
    return p -> sum(f(p) for f in ds)
end
function _compile_waveform(r::PulseRecipe,selected=();base="")
    duration=_compiled_duration(r,base,selected)
    delay=_scalar_accessor(r.delay,_path_join(base,:delay),selected)
    amplitude=_scalar_accessor(r.amplitude,_path_join(base,:amplitude),selected)
    offset=_scalar_accessor(r.offset,_path_join(base,:offset),selected)
    envelope=_compile_envelope(r.envelope,_path_join(base,:envelope),selected,duration)
    modulation=if r.carrier===nothing
        (p,z,t,s) -> z
    else
        frequency=_scalar_accessor(r.carrier.frequency,_path_join(base,"carrier/frequency"),selected)
        phase=_scalar_accessor(r.carrier.phase,_path_join(base,"carrier/phase"),selected)
        reference=r.carrier.reference
        (p,z,t,s) -> begin
            angle=2pi*frequency(p)*(reference===:pulse ? t : s)+phase(p)
            real(z)*cos(angle)-imag(z)*sin(angle)
        end
    end
    return (p,t,sequence_time) -> begin
        l,d=delay(p),duration(p)
        l<=t<=l+d || return offset(p)
        local_time=t==l+d ? d : t-l
        offset(p)+amplitude(p)*modulation(p,envelope(p,local_time),local_time,sequence_time)
    end
end
function _compile_waveform(r::PulseSequenceRecipe,selected=();base="")
    children=map(Tuple(enumerate(r.pulses))) do (i,child)
        prefix=_path_join(base,"pulses/$i")
        (_compile_waveform(child,selected;base=prefix),_compiled_duration(child,prefix,selected),_scalar_accessor(child.delay,_path_join(prefix,:delay),selected))
    end
    duration=_compiled_duration(r,base,selected)
    delay=_scalar_accessor(r.delay,_path_join(base,:delay),selected)
    offset=_scalar_accessor(r.offset,_path_join(base,:offset),selected)
    return (p,t,s) -> begin
        l,d=delay(p),duration(p)
        l<=t<=l+d || return offset(p)
        local_time=t-l; previous=zero(local_time)
        for (i,(f,child_duration,child_delay)) in enumerate(children)
            onset=previous+child_delay(p); stop=onset+child_duration(p)
            if onset<=local_time && (local_time<stop || (i==length(children) && local_time==d))
                return f(p,local_time-previous,local_time)
            end
            previous=stop
        end
        offset(p)
    end
end
