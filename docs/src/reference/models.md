# Components and models

```@docs
Component
make_qubit
make_resonator
make_transmon
make_flux_tunable_transmon
OperatorGenerator
OperatorRegistry
get_operators
get_components
QuantumDevice
add_component!
add_interaction!
add_parameters!
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
