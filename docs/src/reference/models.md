# Components and models

```@docs
Component
make_qubit
make_resonator
make_transmon
make_flux_tunable_transmon
apply_correction
OperatorRegistry
get_operators
get_components
QuantumDevice
add_component!
add_interaction!
Recipe
Leaf
Combine
Input
ResolvedRecipe
draft
graph
add_recipe!
replace_recipe!
NumericalModel
build
states
stateindex
state
stateenergy
```

`get_hamiltonian(device, recipe)` and `get_hamiltonian(model)` return the physical symbolic Hamiltonian;
`make_hamiltonian(model[, free_parameters])`, `make_hamiltonian(prepared)`, and
`make_hamiltonian(model, gate)` return a `QuantumObject` or `QobjEvo`.
With `remove_vacuum=true`, these numerical overloads subtract the model's fixed
reference all-zero occupation energy times identity before `scalar` scaling.
Untracked models use their lowest reference eigenenergy; tracked reductions
missing the vacuum label throw. This is a constant global energy-origin choice,
not an instantaneous ground-state calculation or reference-Hamiltonian inverse.
Stored `model.H` and caches remain unchanged, and evolution/calibration defaults
retain physical phases.

`model.H` is the stored numerical reference. `get_hamiltonian` is symbolic;
use `make_hamiltonian` to copy, compile, scale, or shift a numerical Hamiltonian.

For model free inputs, select symbolic variables or full names from
`model.variables`, for example `make_hamiltonian(model, ["flux_q1", "flux_c"])`.
Each selector is also its runtime dictionary key: string selectors read
`p["flux_q1"]`, and symbolic selectors read `p[flux_q1]`. Mixed selectors are
supported. Each input may appear only once, including across its string and
symbolic forms. With `auto_time_dependence=true`, dictionary entries are functions
of time; with `false`, they are numerical values. Gate pulse-path selectors keep
their existing meaning.
