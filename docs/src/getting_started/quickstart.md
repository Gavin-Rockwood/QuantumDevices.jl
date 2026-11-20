# Build and compose recipes

A `Recipe` defines a model without generating operators or diagonalizing it.
Register the physical components, then select the initial generator dimensions
and the zero-based eigenstates to retain.

```@example quickstart
using QuantumDevices
import Symbolics
using QuantumToolbox

device = QuantumDevice(:three_transmons)
q1 = add_component!(device, make_transmon(:q1, 0.27, 17.0))
q2 = add_component!(device, make_transmon(:q2, 0.29, 19.0))
q3 = add_component!(device, make_transmon(:q3, 0.25, 16.0; dimension=9))
pair = Recipe(device, [q1, q2];
    initial_dims=Dict(q1 => 9, q2 => 9),
    states=Dict(q1 => 0:1, q2 => 0:1),
    interactions=[Interaction(:coupling12,
        0.02 * q1.operators["charge"] * q2.operators["charge"])])
```

The `interactions` keyword copies a collection of `Interaction` objects and binds
their physical component operators to the matching immediate inputs, including
grouped inputs. Every operator must belong to a selected component and name an
exposed operator. Device interactions can be selected by passing their objects or by setting
`interactions=:auto`:

```julia
add_interaction!(device, Interaction(:registered12,
    0.02 * q1.operators["charge"] * q2.operators["charge"]))
automatic_pair = Recipe(device, [q1, q2]; interactions=:auto)
```

Automatic selection includes each registered interaction whose full component
membership is selected and spans at least two immediate inputs. Scalar and
single-component terms are excluded. Terms entirely internal to a grouped input
belong in that upstream recipe, preventing duplicate couplings when composing
recipes. Selection is copied at construction in sorted interaction-name order;
later device registrations do not alter an existing recipe. The default includes
no interactions, and collections do not accept registered interaction names. To add a term
after construction, use `add_interaction!` with that recipe's scoped input
operators, as below.

An explicit input can consume an existing recipe. Its alias selects the grouped
eigenstates to retain, while uncovered components become automatic inputs.
Explicit aliases come first in tensor order, followed by uncovered components
in the order supplied. A grouped input requires the physical component name
between `operators` and the operator name.

```@example quickstart
target = Recipe(device, [q1, q2, q3];
    inputs=(pair=pair,),
    states=Dict(:pair => 0:2, q3 => 0:1))
add_interaction!(target, Interaction(:coupling23,
    0.015 * target.inputs.pair.operators.q2.charge *
        target.inputs.q3.operators.charge))
model = build(device, target)
@assert length(model.spectrum) == 6
@assert all(operator -> operator isa QuantumObject, values(model.operators))
H = make_hamiltonian(model)
```

Only recipes you want to look up by name need registration:

```@example quickstart
add_recipe!(device, :target, target)
rebuilt = build(device, :target; parameters=Dict(q3.variables["EJ"] => 17.0))
```

`initial_dims` applies to automatic component inputs. `states` keys are components
(or their registered symbols) for automatic inputs, and alias symbols for explicit
inputs. Omitted selections keep all states. The default tensor dimension limit
is 1000; `max_dimension=0` disables it. `Leaf`, `Combine`, and `Input` remain
lower-level compatibility interfaces; use `Recipe` for new model definitions.

Hamiltonians use frequency units `E/h`: GHz for time in ns. `make_hamiltonian(model)`
returns a static QuantumToolbox `QuantumObject` in the fixed model construction basis.
Use `make_hamiltonian(model, [parameter])` to obtain a `QobjEvo` whose designated free
parameters read runtime functions from `p[parameter](t)`; set
`auto_time_dependence=false` for scalar runtime values. An empty free-parameter
collection also returns a static `QuantumObject`. Use `make_hamiltonian(model; scalar=2pi)` (or the free-parameter form with the same
keyword) for direct QuantumToolbox evolution. Built-in `evolve` supplies its own
`2π` conversion.

Built models also expose [physical state labels](../user_guide/state_labeling.md)
that follow an interaction ramp while keeping spectra and selections energy-ordered.

Registration captures upstream names, selections, parameter overrides, and local
interactions as an immutable recipe tree. Build an editable trial from a saved
definition:

```julia
trial = draft(device.recipes[:target])
trial.parameters[q3.variables["EJ"]] = 18.0
trial_model = build(device, trial)
add_recipe!(device, :trial, trial)
replace_recipe!(device, :target, trial)
```

A replacement captures a new tree. Existing saved dependencies and model provenance
retain their earlier definitions. Model resolved parameters are immutable snapshot
values. Device components remain device-owned. Arbitrary callables, including
generator closures, can capture mutable external state outside this boundary.
