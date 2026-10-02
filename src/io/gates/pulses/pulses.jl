# Built-in controls are portable recursive records; only custom callables need artifacts.
function _control_record(control::Union{AbstractPulse,AbstractEnvelope,AbstractCarrier}, directory, index)
    throw(ArgumentError("No persistence rule for $(typeof(control))"))
end
_control_record(::Constant, directory, index) = Dict("type" => "constant")
_control_record(::SineSquared, directory, index) = Dict("type" => "sine_squared")
function _control_fields(control, directory, index; skip=())
    return Dict(String(key) => _encode_control_value(getfield(control, key), directory, "$index-$key")
        for key in fieldnames(typeof(control)) if !(key in skip))
end
function _typed_control_record(tag, control, directory, index)
    return Dict("type" => tag, "fields" => _control_fields(control, directory, index))
end
_control_record(x::Pulse, directory, index) = _typed_control_record("pulse", x, directory, index)
_control_record(x::Gaussian, directory, index) = _typed_control_record("gaussian", x, directory, index)
_control_record(x::RampedFlattop, directory, index) = _typed_control_record("ramped_flattop", x, directory, index)
_control_record(x::SineCarrier, directory, index) = _typed_control_record("sine_carrier", x, directory, index)
function _custom_control_record(tag, control, directory, index)
    mkpath(directory)
    filename = "$index.jld2"
    JLD2.save_object(_bundle_child(directory, filename), control.callable)
    return Dict("type" => tag, "fields" => _control_fields(control, directory, index; skip=(:callable,)),
        "artifact" => filename)
end
_control_record(x::Envelope, directory, index) = _custom_control_record("custom_envelope", x, directory, index)
_control_record(x::Carrier, directory, index) = _custom_control_record("custom_carrier", x, directory, index)

_encode_control_value(x, directory, index) = _encode_value(x)
_encode_control_value(x::Union{AbstractPulse,AbstractEnvelope,AbstractCarrier}, directory, index) =
    Dict("kind" => "control", "record" => _control_record(x, directory, index))
_encode_control_value(x::NamedTuple, directory, index) = Dict("kind" => "named_tuple",
    "fields" => Dict(String(k) => _encode_control_value(v, directory, "$index-$k") for (k,v) in pairs(x)))
function _decode_control_value(data, directory)
    if !(data isa Union{Nothing,Number,AbstractString})
        kind = data["kind"]
        kind == "control" && return _restore_control(data["record"], directory)
        kind == "named_tuple" && return (; (Symbol(k) => _decode_control_value(v, directory)
            for (k,v) in pairs(data["fields"]))...)
    end
    return _decode_value(data)
end
function _restore_control_fields(data, directory)
    return (; (Symbol(k) => _decode_control_value(v, directory) for (k,v) in pairs(data["fields"]))...)
end
_restore_control(data, directory) = _restore_control(Val(Symbol(data["type"])), data, directory)
_restore_control(::Val{T}, data, directory) where {T} =
    throw(ArgumentError("Unknown or legacy pulse record: $T; recreate it with Pulse"))
_restore_control(::Val{:constant}, data, directory) = Constant()
_restore_control(::Val{:sine_squared}, data, directory) = SineSquared()
function _restore_control(::Val{:pulse}, data, directory)
    f = _restore_control_fields(data, directory)
    return Pulse(f.envelope; carrier=f.carrier, amplitude=f.amplitude, offset=f.offset,
        duration=f.duration, delay=f.delay)
end
function _restore_control(::Val{:gaussian}, data, directory)
    f = _restore_control_fields(data, directory)
    return Gaussian(f.sigma; center=f.center)
end
function _restore_control(::Val{:ramped_flattop}, data, directory)
    f = _restore_control_fields(data, directory)
    return RampedFlattop(f.rise_time, f.fall_time, f.ramp_up, f.ramp_down, f.split_up, f.split_down)
end
function _restore_control(::Val{:sine_carrier}, data, directory)
    f = _restore_control_fields(data, directory)
    return SineCarrier(f.frequency; phase=f.phase, reference=f.reference)
end
function _restore_control(::Val{:custom_envelope}, data, directory)
    f = _restore_control_fields(data, directory)
    callable = JLD2.load_object(_bundle_child(directory, data["artifact"]))
    return Envelope(callable, f.parameters)
end
function _restore_control(::Val{:custom_carrier}, data, directory)
    f = _restore_control_fields(data, directory)
    callable = JLD2.load_object(_bundle_child(directory, data["artifact"]))
    return Carrier(callable, f.parameters; reference=f.reference)
end
