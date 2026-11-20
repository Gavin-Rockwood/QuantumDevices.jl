using JLD2
# New identifiers describe normalized payloads; v1 identifiers retain physical semantics.
const _PULSE_IMPLEMENTATIONS=(pulse=Pulse,pulse_sequence=PulseSequence,constant_frequency=ConstantFrequency,
    normalized_envelope=Envelope,normalized_constant=Constant,normalized_gaussian=Gaussian,
    normalized_gaussian_zero=GaussianZero,normalized_sine_squared=SineSquared,normalized_cosine=Cosine,
    normalized_blackman=Blackman,normalized_sech=Sech,normalized_bump=Bump,normalized_drag=DRAG,
    normalized_slepian=Slepian,normalized_gaussian_square=GaussianSquare,normalized_erf_square=ErfSquare,
    normalized_flattop=RampedFlattop,normalized_interpolation=InterpolatedEnvelope,
    normalized_region=Region,normalized_custom=CustomShape,legacy_physical=LegacyEnvelope,
    envelope_sequence=EnvelopeSequence,segment=Segment)
const _V1_SHAPES=(constant=:Constant,gaussian=:Gaussian,gaussian_zero=:GaussianZero,sine_squared=:SineSquared,
    cosine=:Cosine,blackman=:Blackman,sech=:Sech,bump=:Bump,drag=:DRAG,slepian=:Slepian,
    gaussian_square=:GaussianSquare,erf_square=:ErfSquare,ramped_flattop=:RampedFlattop,
    ramp=:Ramp,interpolated_envelope=:InterpolatedEnvelope,custom_envelope=:Envelope)
function _pulse_storage_fields(x)
    x isa LegacyEnvelope && return (:implementation,:fields,:duration)
    x isa PulseSequence && return (:pulses,:offset,:delay)
    x isa Slepian && return (:time_bandwidth,:samples)
    return fieldnames(typeof(x))
end
function _pulse_record(x::_PulseObject)
    tag=findfirst(T->x isa T,_PULSE_IMPLEMENTATIONS)
    tag===nothing && throw(ArgumentError("No persistence implementation for $(typeof(x))"))
    return (kind=:control,implementation=tag,version=2,
        fields=(;(k=>_pulse_record(getproperty(x,k)) for k in _pulse_storage_fields(x))...))
end
function _restore_pulse_fields(id,version,f::NamedTuple)
    version in (1,2) || throw(ArgumentError("Unsupported pulse record version"))
    if version==1 && haskey(_V1_SHAPES,id)
        kind=_V1_SHAPES[id]
        T=_LEGACY_SHAPE_TYPES[kind]
        expected=filter(k->!(k in (:interpolant,:weights)),fieldnames(T))
        Set(keys(f))==Set(expected) || throw(ArgumentError("Legacy pulse field schema mismatch"))
        return LegacyEnvelope(kind,NamedTuple{Tuple(expected)}(Tuple(getproperty(f,k) for k in expected)),nothing)
    end
    version==1 && !(id in (:pulse,:pulse_sequence,:constant_frequency,:envelope_sequence,:segment)) &&
        throw(ArgumentError("Normalized shape and adapter identifiers require pulse record version 2"))
    haskey(_PULSE_IMPLEMENTATIONS,id) || throw(ArgumentError("Unknown pulse implementation $id"))
    T=_PULSE_IMPLEMENTATIONS[id]
    expected=T===Slepian ? (:time_bandwidth,:samples) : fieldnames(T)
    T===LegacyEnvelope && (expected=(:implementation,:fields,:duration))
    version==1 && T===Pulse && (expected=(:envelope,:carrier,:amplitude,:offset,:duration,:delay))
    version==1 && T===Segment && (expected=(:shape,:duration))
    T===PulseSequence && (expected=(:pulses,:offset,:delay))
    Set(keys(f))==Set(expected) || throw(ArgumentError("Pulse field schema mismatch for $id"))
    T===Pulse && return Pulse(f.envelope,f.carrier,f.amplitude,f.offset,version==1 ? f.duration : nothing,f.delay)
    T===PulseSequence && return PulseSequence(f.pulses...;offset=f.offset,delay=f.delay)
    T===EnvelopeSequence && return EnvelopeSequence(f.segments...)
    T===Segment && return Segment(f.shape;duration=f.duration,scale=get(f,:scale,1))
    T===Envelope && return _raw_envelope(f.shape,f.duration,f.rise,f.fall)
    T===Slepian && return Slepian(Val(:raw),f.time_bandwidth,f.samples)
    T in (Constant,SineSquared,Cosine,Blackman) && return T(Val(:raw))
    if T in (Gaussian,GaussianZero,DRAG,Sech,Bump,GaussianSquare,ErfSquare,RampedFlattop,InterpolatedEnvelope)
        T===InterpolatedEnvelope && !(f.times isa Tuple && f.values isa Tuple) &&
            throw(ArgumentError("Normalized interpolation samples must be immutable tuples"))
        shape=T(Val(:raw),(getproperty(f,k) for k in expected)...)
        _validate_shape(shape)
        return shape
    end
    if T===Region
        shape=Region(f.shape,_normalized_interval(f.interval))
        _validate_shape(shape)
        return shape
    end
    shape=T((getproperty(f,k) for k in expected)...)
    shape isa _ShapePayload && _validate_shape(shape)
    return shape
end
_pulse_record(x::NamedTuple) =
    (kind = :named_tuple, fields = (; (k=>_pulse_record(v) for (k, v) in pairs(x))...))
_pulse_record(x::Tuple) = (kind = :tuple, values = map(_pulse_record, x))
_pulse_record(x::AbstractVector) = (kind = :vector, values = map(_pulse_record, x))
_pulse_record(x) = (kind = :value, value = x)
function _restore_pulse_record(r)
    r.kind === :value && return r.value
    r.kind === :tuple && return map(_restore_pulse_record, r.values)
    r.kind === :vector && return map(_restore_pulse_record, r.values)
    f=(; (k=>_restore_pulse_record(v) for (k, v) in pairs(r.fields))...)
    r.kind === :named_tuple && return f
    r.kind === :control || throw(ArgumentError("Unknown pulse record"))
    return _restore_pulse_fields(r.implementation,get(r,:version,1),f)
end
"""
    save_pulse(filename, control)

Save a versioned recursive pulse/envelope record with JLD2. Builtins store stable
implementation identifiers and constructor parameters, never interpolation caches.
Custom `Envelope` callables use JLD2's function serialization; their definitions
and dependencies must be available when loading, and anonymous functions may not
be portable between sessions. Only load files from trusted sources.
"""
function save_pulse(filename::AbstractString, control::_PulseObject)
    JLD2.jldsave(filename; pulse_schema = 1, control = _pulse_record(control))
    return filename
end
"""Load a pulse record and rebuild each object through its validated constructor."""
function load_pulse(filename::AbstractString)
    JLD2.jldopen(filename, "r") do file
        file["pulse_schema"] == 1 || throw(ArgumentError("Unsupported pulse file schema"))
        return _restore_pulse_record(file["control"])
    end
end
