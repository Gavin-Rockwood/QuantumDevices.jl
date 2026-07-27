function _device_names(registry)
    join(sort!(collect(keys(registry))), ", ")
end

_device_names_or_none(registry) =
    isempty(registry) ? "none" : _device_names(registry)

function Base.show(io::IO, device::QuantumDevice)
    print(io, "QuantumDevice(:", device.name,
        "; components=[", _device_names(device.components),
        "], interactions=[", _device_names(device.interactions),
        "], models=[", _device_names(device.modelspecs),
        "], gates=[", _device_names(device.gatespecs), "])")
end

function Base.show(io::IO, ::MIME"text/plain", device::QuantumDevice)
    println(io, "QuantumDevice :", device.name)
    println(io, "  Components: ", _device_names_or_none(device.components))
    println(io, "  Interactions: ", _device_names_or_none(device.interactions))
    println(io, "  Models: ", _device_names_or_none(device.modelspecs))
    print(io, "  Gates: ", _device_names_or_none(device.gatespecs))
end
