# Parameter discovery and updates

`parameters(gate)` returns `name => (path, value)`. Gate scalars retain their names;
pulse parameters use names such as `"drive/amplitude"`, `"drive/carrier/frequency"`,
and `"drive/envelope/sigma"`, including nested ramp and custom shape parameters. `calibration_values` is an
alias; neither function filters the discovery map to optimizable numbers.

```@example parameters
using QuantumDevices
pulse = Pulse(Gaussian(0.1); amplitude=0.2, duration=1.0)
gate = DeviceGate((; drive=pulse), param(:drive) * op(:q_x), 1.0)
path, original = parameters(gate)["drive/amplitude"]
updated = setpath(gate, path, 0.4)
@assert getpath(updated, path) == 0.4
@assert getpath(gate, path) == original
@assert haspath(gate, "drive/envelope/sigma")
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

Reconstruction follows the changed inputs:

| Update | Recomputed | Reused |
|---|---|---|
| Gate parameters | Updated gate and its propagator when evolved | Idle Hamiltonian, component bases, dressed states |
| Interaction strengths | Interaction coefficients, idle Hamiltonian, dressed-state tracking | Local eigensystems, uncoupled preparation, projected operator terms |
| Interaction expression | Interaction operator terms and Hamiltonian, dressed-state tracking | Component bases and uncoupled preparation |
| Component physics, such as transmon `EC` or `EJ` | Affected component basis/operators, dependent projections and Hamiltonian terms, dressed-state tracking | Unchanged components and their preparation |
| Dimensions, truncations, or dimension limit | Complete model reconstruction | Existing component objects when applicable |

A transmon physical update reconstructs its parent charge-basis eigensystem and
transforms its operators before projection to the retained space. Products such
as `charge^2` are recomputed in the parent space before projection. Other
components do not get rediagonalized for this update.

Derived properties cannot be replaced directly. Component names and order must
survive a model component update. Supported updates do not change the original
model; unchanged branches and read-only preparation can be shared. Use these
updates rather than mutating cached matrices or physical input containers.

The reconstruction cache is derived and is rebuilt when loading a saved model;
the bundle schema does not change.

Custom pulse types can expose selectable values through `parameters(pulse)`;
see the [extension guide](../development/extensions.md).
