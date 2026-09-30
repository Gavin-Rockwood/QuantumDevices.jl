# Gates and evolution

A `DeviceGate` contains parameters, an additional Hamiltonian, and duration.
Parameters must be finite numbers or pulses; raw functions are rejected. Use
`GenericPulseFunction(f, parameters)` for user-defined `f(p, t)` controls.

```@example gates
using QuantumDevices, QuantumToolbox, SciMLBase, LinearAlgebra
q = make_qubit("q", 0.0)
model = make_model([q], val(0), (;))
gate = DeviceGate((; drive=constant_pulse(π / 2)),
                  param(:drive) * op(:q_x), 1.0)
U = gate_unitary(model, gate; abstol=1e-9, reltol=1e-9)
@assert unitary_infidelity(-im * sigmax(), U) < 1e-7
unitary_infidelity(-im * sigmax(), U)
```

The ideal zero-idle example isolates an X rotation. In a nonzero idle model, lab
frame detuning and drift remain present unless explicitly modeled or compensated.
There is no automatic rotating frame or rotating-wave approximation.

For state evolution, use the same Hamiltonian with QuantumToolbox:

```@example gates
H = numerical(model, gate)
result = sesolve(H, basis(2, 0), range(0, 1; length=21); progress_bar=false)
@assert SciMLBase.successful_retcode(result.retcode)
@assert abs2(dot(basis(2, 1), result.states[end])) > 0.999
length(result.states)
```

## Overrides and endpoint checks

Gate parameters override matching model parameters. A parameter that changes an
idle coefficient must equal the idle value **exactly** at both endpoints; this is
an equality check, not a tolerance-based continuity test. An additional drive
parameter absent from the idle model has no idle-matching restriction.

Duration must be finite and nonnegative; `gate_unitary` requires positive duration.
A zero additional Hamiltonian represents a gate that changes idle coefficients
only. `numerical` returns a static quantum operator when all coefficients are
static, and a `QobjEvo` when time dependence is present.

Named gates are stored in `model.gates[:name]`. A gate object is immutable;
use [parameter paths](parameters.md) to construct updated versions.
