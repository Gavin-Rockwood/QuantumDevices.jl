function pulse_duration(p::Union{Pulse,AbstractPulseRecipe},resolved::NamedTuple)
    r=p isa Pulse ? p.recipe : p
    r isa PulseRecipe && return resolved.duration
    return foldl(_checked_time_add,(_checked_time_add(child.delay,pulse_duration(r.pulses[i],child)) for (i,child) in enumerate(resolved.pulses)))
end
# Timing consumes the same final physical fields as waveform evaluation.
function _resolved_shape_stops(s,p,d)
    s isa InterpolatedEnvelope && return collect(p.times[2:end-1])
    s isa Union{GaussianSquare,ErfSquare} && return iszero(p.ramp_time) ? Real[] : [p.ramp_time,d-p.ramp_time]
    if s isa Region
        a,b=s.interval; source_duration=d/(b-a)
        return [t-a*source_duration for t in _resolved_shape_stops(s.shape,p,source_duration) if a*source_duration<t<b*source_duration]
    elseif s isa RampedFlattop
        stops=Real[p.rise_time,d-p.fall_time]
        if s.ramp_up isa NamedTuple && p.rise_time>0
            append!(stops,_resolved_shape_stops(s.ramp_up.shape,p.ramp_up,p.rise_time))
        end
        if s.ramp_down isa NamedTuple && p.fall_time>0
            append!(stops,d-p.fall_time .+ _resolved_shape_stops(s.ramp_down.shape,p.ramp_down,p.fall_time))
        end
        return sort!(unique(filter(t -> 0<t<d,stops)))
    end
    return [u*d for u in _shape_tstops(s)]
end
_resolved_envelope_stops(e::Envelope,p,d)=_resolved_shape_stops(e.shape,p,d)
function _resolved_envelope_stops(e::EnvelopeSequence,p,d)
    stops=Real[]; start=zero(d)
    for (s,v) in zip(e.segments,p.segments)
        append!(stops,start .+ _resolved_envelope_stops(s.shape,v.shape,v.duration))
        start+=v.duration; push!(stops,start)
    end
    return stops
end
function pulse_tstops(p::Union{Pulse,AbstractPulseRecipe},resolved::NamedTuple)
    r=p isa Pulse ? p.recipe : p
    result=Real[resolved.delay,resolved.delay+pulse_duration(r,resolved)]
    if r isa PulseRecipe
        append!(result,resolved.delay .+ _resolved_envelope_stops(r.envelope,resolved.envelope,resolved.duration))
    else
        start=resolved.delay
        for (child,v) in zip(r.pulses,resolved.pulses)
            append!(result,start .+ pulse_tstops(child,v))
            start+=v.delay+pulse_duration(child,v)
        end
    end
    return sort!(unique(result))
end
_resolved_envelope_stops(e::LegacyEnvelope,p,d)=_validated_envelope_tstops(e,d)

_resolved_envelope_stops(e::AbstractEnvelope,p,d)=_validated_envelope_tstops(e,d)
