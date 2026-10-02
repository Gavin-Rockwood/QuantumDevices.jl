"""
    gate_infidelity(model::DeviceModel, gate::DeviceGate, target)

Evaluate [`unitary_infidelity`](@ref) between `target` and
[`get_unitary`](@ref). The target acts on the full retained model space.
This is the default calibration objective, not a computational-subspace
leakage metric or open-system fidelity.
"""
gate_infidelity(model::DeviceModel, gate::DeviceGate, target) =
    unitary_infidelity(target, get_unitary(model, gate))
