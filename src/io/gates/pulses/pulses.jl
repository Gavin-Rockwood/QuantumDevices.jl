# Extend these two methods for additional pulse subtypes and their storage tags.
_pulse_record(pulse::AbstractPulse, directory, index) =
    throw(ArgumentError("No persistence rule for $(typeof(pulse))"))
_restore_pulse(::Val{T}, data, directory) where {T} =
    throw(ArgumentError("Unknown pulse storage tag: $T"))
_restore_pulse(data, directory) = _restore_pulse(Val(Symbol(data["type"])), data, directory)

include("InternalPulseFunction.jl")
include("GenericPulseFunction.jl")
