# Parameter discovery and updates

`parameters(gate)` returns `name => (path, value)`. Gate scalars retain their names;
pulse parameters use names such as `"drive/amplitude"`. `calibration_values` is an
alias; neither function filters the discovery map to optimizable numbers.

```@example parameters
using QuantumDevices
pulse = gaussian_pulse(0.2, 0.1)
gate = DeviceGate((; drive=pulse), param(:drive) * op(:q_x), 1.0)
path, original = parameters(gate)["drive/amplitude"]
updated = setpath(gate, path, 0.4)
@assert getpath(updated, path) == 0.4
@assert getpath(gate, path) == original
@assert haspath(gate, "drive/sigma")
(path, original)
```

Paths use `/` separators, not Julia expressions. They traverse fields, named tuples,
dictionaries, arrays, and tuples. Indices are one-based. A parameter lookup can omit
`parameters/` unless a structure field shadows it; explicit paths always resolve
that ambiguity. A dictionary containing both `:a` and `"a"` is ambiguous for `"a"`.

`setpath` copies the edited path and rebuilds immutable structs; it is not a deep
copy of every unchanged branch. `setpath!` requires a mutable root. A model or gate
root is immutable, although `model.gates` is a mutable dictionary.

## Physical reconstruction

Updating built-in component parameters reconstructs its basis and derived
parameters. Updating a model's components, coupling parameters, interactions,
truncations, or dimension limit reconstructs its Hamiltonian and dressed states.
Changing its gates preserves the idle physics. Derived properties cannot be
replaced directly. Component names and order must survive a model component update.

Custom pulse types can expose selectable values through `parameters(pulse)`;
see the [extension guide](../development/extensions.md).
