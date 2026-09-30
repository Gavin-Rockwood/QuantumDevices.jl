# Your first model and gate

This example constructs two ideal qubits, evaluates a shaped drive, and evolves a
state. The chosen coefficients and time are dimensionless consistent units.

```@example quickstart
using QuantumDevices, QuantumToolbox, LinearAlgebra
q1 = make_qubit("q1", 1.0)
q2 = make_qubit("q2", 1.4)
model = make_model([q1, q2], param(:g) * op(:q1_x) * op(:q2_x), (; g=0.02))
@assert size(model.H) == (4, 4)
size(model.H)
```

The model operator names are prefixed with component names. The interaction
couples their Pauli X operators. The component order is `q1 ⊗ q2`.

```@example quickstart
pulse = ramped_flattop_pulse(0.1, 0.2)
gate = DeviceGate((; drive=pulse), param(:drive) * op(:q1_x), 1.0)
H = numerical(model, gate)
@assert H(0.5) ≈ model.H + 0.1 * model.operators.q1_x
size(H(0.5))
```

The drive adds to the idle Hamiltonian. Its frequency is not automatically chosen
to be resonant: this is an example of building a control, not an X-gate calibration.

```@example quickstart
ψ0 = tensor(basis(2, 0), basis(2, 0))
solution = sesolve(H, ψ0, range(0, 1; length=21); progress_bar=false)
@assert isapprox(norm(solution.states[end]), 1; atol=1e-6)
length(solution.states)
```

Store a named gate and save a bundle to a new destination:

```@example quickstart
model.gates[:drive] = gate
mktempdir() do directory
    path = save(joinpath(directory, "device"), model)
    restored = load(path)
    @assert restored.H ≈ model.H
    @assert numerical(restored, restored.gates[:drive])(0.5) ≈ H(0.5)
    println("Model and gate restored successfully")
end
```

Next: [components](../user_guide/components.md), [models](../user_guide/models.md),
and [pulses](../user_guide/pulses.md).
