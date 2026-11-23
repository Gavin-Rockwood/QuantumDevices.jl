# Register only fully validated immutable definitions; failed replacement is atomic.
"""Register an owned static map after checking every output belongs to the device. Names must be new."""
function add_control_map!(device::QuantumDevice, map::AbstractControlMap)
    haskey(device.control_maps,map.name) && throw(ArgumentError("control map $(map.name) already registered"))
    _validate_control_map_outputs(map,Iterators.flatten(values(c.variables) for c in values(device.components)))
    device.control_maps[map.name] = _copy_control_map(map)
    return device.control_maps[map.name]
end
"""Atomically replace an existing registry map. Captured recipe and gate maps retain their values."""
function replace_control_map!(device::QuantumDevice, map::AbstractControlMap)
    haskey(device.control_maps,map.name) || throw(ArgumentError("unknown control map $(map.name)"))
    _validate_control_map_outputs(map,Iterators.flatten(values(c.variables) for c in values(device.components)))
    device.control_maps[map.name] = _copy_control_map(map)
    return device.control_maps[map.name]
end
