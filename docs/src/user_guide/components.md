# Components

A component declares symbolic scalar variables, their idle values, an owner-aware
Hamiltonian, and numerical operator factories. Construction does not diagonalize
or generate a numerical operator. Recipes select dimensions and retained states;
`build` prepares the numerical model.

```@example components
using QuantumDevices
q = make_qubit(:q; frequency=1.0)
r = make_resonator(:r; frequency=2.5, dimension=4)
t = make_transmon(:t, 0.2, 5.0; dimension=9)
ft = make_flux_tunable_transmon(:ft, 0.2, 9.0; dimension=9, flux=0.2)
@assert t.parameters["EC"] == 0.2
@assert Set(keys(ft.variables)) == Set(["EC","EJmax","d","ng","flux"])
ft.hamiltonian
```

| Component | Local variable aliases | Local operators | Numerical parent basis |
|:--|:--|:--|:--|
| Qubit | `ν` | `x`, `y`, `z` | Two-state Pauli basis |
| Resonator | `ν` | `a`, `n` | Truncated Fock basis |
| Transmon | `EC`, `EJ`, `ng` | `jump`, `charge` | Finite charge basis |
| Flux-tunable transmon | `EC`, `EJmax`, `d`, `ng`, `flux` | `jump`, `charge` | Finite charge basis |

`component.variables` maps local strings to atomic Symbolics variables.
`component.parameters` maps the same strings to finite concrete idle numbers.
`component.operators` maps local strings to symbolic operators, and
`component.operator_generators` maps the same strings to numerical factories.
Aliases are local to each component. By default symbolic names use its name as
a suffix, such as `EC_q1`. All builtin constructors accept `subscript=missing`
for that default, `subscript=nothing` for no suffix, or a nonempty Symbol/string
for a custom suffix. The display `label` is independent. Equal symbolic names
share a parameter and must have consistent idle values.
Use the actual variables as recipe override keys and control targets. Individual
variables carry no idle-value field.

## Dimensions and convergence

Resonator dimensions must be positive integers. Transmon dimensions must be
positive odd integers: `dimension = 2ncut + 1`, covering charge states
`-ncut:ncut`. Increase the parent cutoff and retained-state count independently
to check convergence. `initial_dims` in a recipe can override factory dimensions.
All transmon Hamiltonian inputs use frequency units `E/h`, and `flux` is `Φext/Φ0`.

`build(device, recipe; parameters=Dict(t.variables["EJ"]=>6.0))` recomputes the
physical basis at that operating point. `make_hamiltonian(model,[ft.variables["flux"]])`
varies only the scalar coefficients in the existing fixed basis.

## Custom components

Use `component_parameters(name; alias=value, ...)` for the variable and idle
maps, and `component_operators(name; local_name=OperatorGenerator(...))` for
the symbolic operator and factory maps. Pass all four maps to
`Component(name; variables, parameters, operators, operator_generators, hamiltonian)`.
The factory callback accepts a requested dimension and any keyword options.
See [extensions](../development/extensions.md) for a complete example.

## Device operator ownership and interaction membership

Registering a component freezes its symbolic operator declarations and generator keys.
Generator values may still be replaced under an existing key.
The component keeps its identity, so both a component passed to `add_component!`
and the returned component can be passed to `Recipe`. Draft components remain
editable. Constructor input dictionaries are copied; changing them cannot change
registered declarations. Parameter defaults remain editable and are validated
when the device defaults are assembled. Generator callables retain their identity;
mutable state captured by a callable remains the caller's responsibility.

```julia
device = QuantumDevice(:pair)
q1 = add_component!(device, make_qubit(:q1))
q2 = add_component!(device, make_qubit(:q2))
device.operators[:q1]["x"]             # the physical symbolic operator
device.operators[q1.operators["x"]]    # :q1
coupling = add_interaction!(device,
    Interaction(:coupling, 0.01 * q1.operators["x"] * q2.operators["x"]))
device.interactions_for[:q1]           # (:coupling,)
device.interaction_components[:coupling] # (:q1, :q2)
get_operators(coupling, device)         # unique physical operators
get_components(coupling, device)       # [:q1, :q2]
```

The registry also accepts an `OperatorRef`, a scoped recipe operator, or its bound
reference, resolving each to its individual physical operator. Composite operator
expressions are rejected as registry keys. Missing registry or membership lookups
raise `ArgumentError`. Inspection helpers preserve first operator occurrence order;
without a device they inspect physical identities, and with a device they also
validate registration.

Membership values are immutable tuples of names, resolved through
`device.interactions` or `device.components`. Each interaction is one hyperedge,
including interactions involving three or more components. Scalar interactions
have an empty component tuple, and a single-component interaction has one owner.
Recipe-local interactions never enter these device indexes. Registered interactions
capture their operator-relevant definition; numerical defaults remain editable.
Register components and interactions through `add_component!` and
`add_interaction!` so validation and all indexes update together.

`QuantumDevice(saved_context)` rebuilds independent physical definitions and both indexes
from a saved gate or calibration context.

`device.interactions_for[q1]` and `device.operators[q1]` also accept the registered
component object, equivalent to indexing by `:q1`. A separate same-name component
is rejected; use the object registered on that device (including after restoring
a saved context). Both indexes have compact summaries and grouped notebook displays;
interaction membership displays show each interaction's full set of component owners.
