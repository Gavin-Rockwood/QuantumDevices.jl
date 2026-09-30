# Symbolic Hamiltonians

Use `op` for named operators, `param` for scalar parameters, and `val` for numeric
constants. Addition, ordered multiplication, powers, and supported scalar functions
form expressions. `call(f, args...)` can express a scalar coefficient explicitly.

```@example symbolics
using QuantumDevices, QuantumToolbox
expression = param(:a) * op(:x) + param(:b) * op(:z)
operators = (; x=sigmax(), z=sigmaz())
H = numerical(expression, operators, (; a=0.2, b=0.5))
@assert H ≈ 0.2 * sigmax() + 0.5 * sigmaz()
size(H)
```

Supplying all scalar parameters produces a static result. Missing parameters or
function-valued parameters produce a `QobjEvo`. Dynamic scalar functions receive
physical time as their single argument.

```@example symbolics
dynamic = numerical(expression, operators, (; b=0.5))
@assert dynamic((; a=0.2), 0.0) ≈ H
timed = numerical(expression, operators, (; a=t -> 0.2cos(t), b=0.5))
@assert timed(0.0) ≈ H
size(timed(0.0))
```

At the low-level expression interface, unresolved parameters can be supplied when
calling the evolution object. Models require static idle parameters; gates require
all referenced parameters and use pulse objects for controls.

## Names and products

A local component Hamiltonian uses names such as `op(:charge)` and `param(:EC)`.
The model promotes these to `op(:t_charge)` and `param(:t_EC)` for component `"t"`.
Interactions must use these promoted names. Coupling parameters retain their names.

Operator words preserve order. General operator-valued function calls do not form
part of the sum/product compiler. Arbitrary scalar callables may evaluate but
cannot necessarily be saved: symbolic persistence uses an operation allowlist.
See the [symbolic API](../reference/symbolics.md).
