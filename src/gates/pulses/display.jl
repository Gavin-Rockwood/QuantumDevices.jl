Base.show(io::IO,e::LegacyEnvelope)=print(io,"LegacyEnvelope{",e.implementation,"}(duration=",e.duration,")")
function Base.show(io::IO,e::Envelope)
    name=e.shape isa Region ? "Region{$(nameof(typeof(e.shape.shape)))}" : string(nameof(typeof(e.shape)))
    print(io,"Envelope{",name,"}(duration=",e.duration)
    for (key,value) in pairs(_physical_parameters(e.shape,e.duration))
        value isa Number && print(io,", ",key,"=",value)
    end
    return print(io,")")
end
Base.show(io::IO, p::Pulse) =
    print(io, "Pulse(", p.envelope, ", duration=", p.duration, ", delay=", p.delay, ")")
Base.show(io::IO, p::PulseSequence) = print(
    io,
    "PulseSequence(",
    length(p.pulses),
    " children, duration=",
    p.duration,
    ", delay=",
    p.delay,
    ")",
)
