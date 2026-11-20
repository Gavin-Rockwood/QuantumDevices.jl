function _pulse_value(p::Pulse, t, sequence_time)
    p.delay<=t<=p.delay+p.duration || return p.offset
    local_time=t==p.delay+p.duration ? p.duration : t-p.delay
    return _require_finite(
        "Pulse value",
        p.offset+p.amplitude*_modulate(
            p.carrier,
            p.envelope(local_time, p.duration),
            local_time,
            sequence_time,
        ),
    )
end
function _pulse_value(p::PulseSequence, t, sequence_time)
    p.delay<=t<=p.delay+p.duration || return p.offset
    local_time=t-p.delay
    previous=zero(local_time)
    for (i, child) in enumerate(p.pulses)
        onset=previous+child.delay
        stop=onset+child.duration
        if onset<=local_time && (local_time<stop || (i==length(p.pulses) && local_time==p.duration))
            child_time=local_time==stop ? child.delay+child.duration : local_time-previous
            return _pulse_value(child, child_time, local_time)
        end
        previous=stop
    end
    return p.offset
end
function (p::AbstractPulse)(t::Real)
    _require_resolved_pulse(p)
    _require_finite_real("time", t)
    return _pulse_value(p, t, t)
end
"""
    pulse_tstops(control)

Return sorted unique physical control boundaries, joins and interpolation knots,
including delay and nested sequence placement. Supply these times to an integrator.
"""
pulse_tstops(p::Pulse) = sort!(
    unique(
        vcat([p.delay, p.delay+p.duration], p.delay .+ _validated_envelope_tstops(p.envelope, p.duration)),
    ),
)
function pulse_tstops(p::PulseSequence)
    result=Real[p.delay, p.delay+p.duration]
    start=p.delay
    for child in p.pulses
        append!(result, start .+ pulse_tstops(child))
        start+=child.delay+child.duration
    end
    return sort!(unique(result))
end

# Offsets alone may name a captured model reference. Never consult variable
# metadata, pulse targets, or interaction defaults when resolving a baseline.
function _validate_pulse_offset(offset)
    if offset isa Num
        is_parameter(offset) || throw(ArgumentError("pulse offset must be numeric or an atomic Symbolics reference"))
    else
        _require_finite("offset", offset)
    end
    return offset
end
function _reference_offset(offset, reference)
    offset isa Num || return offset
    haskey(reference, offset) || throw(ArgumentError("pulse offset reference $offset is missing from model.parameters"))
    return reference[offset]
end
function _resolve_pulse(p::Pulse, reference)
    return Pulse(p.envelope, p.carrier, p.amplitude,
        _reference_offset(p.offset, reference), p.duration, p.delay)
end
_resolve_pulse(p::PulseSequence, reference) = PulseSequence(
    (_resolve_pulse(child, reference) for child in p.pulses)...;
    offset=_reference_offset(p.offset, reference), delay=p.delay)
function _require_resolved_pulse(p::AbstractPulse)
    p.offset isa Num && throw(ArgumentError("pulse offset $(p.offset) requires a model reference; compile with make_hamiltonian(model, gate) or build(device, gate)"))
    p isa PulseSequence && foreach(_require_resolved_pulse, p.pulses)
    return nothing
end
