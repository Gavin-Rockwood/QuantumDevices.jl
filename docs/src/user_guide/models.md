# Coupled models and truncation

`make_model` combines components in the supplied tensor order and compiles their
local Hamiltonians and interaction. Component names must be unique, and coupling
parameters must not collide with promoted component parameters.

```@example models
using QuantumDevices, QuantumToolbox
q = make_qubit("q", 1.0)
t = make_transmon("t", 0.2, 5.0, 9)
interaction = param(:g) * op(:q_x) * op(:t_charge)
model = make_model([q, t], interaction, (; g=0.01);
                   truncation_dimensions=Dict(t => 3), max_dimension=20)
@assert size(model.H) == (6, 6)
@assert model.parameters.t_EC == 0.2
size(model.operators.t_charge)
```

`truncation_dimensions` keys are the component objects, not names. Each retained
dimension must be between one and its parent dimension. The product is bounded by
`max_dimension`; this check prevents accidental construction of an oversized space,
but is not a memory-use prediction for all subsequent solvers.

## Model properties

| Property | Meaning |
|:--|:--|
| `components` | Original components in tensor order |
| `parameters` | Merged static component and coupling values |
| `operators` | Individually projected and embedded named operators |
| `hamiltonian` | Complete symbolic Hamiltonian |
| `H` | Compiled sparse idle Hamiltonian |
| `states` | Dressed states labeled by bare product-level tuples |
| `others` | Dressed energies with the same labels |
| `confidence` | Continuity diagnostic for the state labels |
| `gates` | Mutable dictionary of named gate objects |

A complete operator word is multiplied in the parent basis **before** projection.
Consequently, `numerical` of `op(:t_charge)^2` can differ from squaring
`model.operators.t_charge`. Use symbolic gate/model expressions when the distinction
matters. See [projection](../explanations/projection.md).

Use `numerical(model, expression)` to evaluate a symbolic expression with the
model's parameters and projection rules:

```@example models
H = numerical(model, model.hamiltonian)
@assert H ≈ model.H
charge_squared = numerical(model, op(:t_charge)^2)
size(charge_squared)
```

## Changing inputs

```@example models
updated = setpath(model, "coupling_parameters/g", 0.02)
@assert model.coupling_parameters.g == 0.01
@assert updated.coupling_parameters.g == 0.02
@assert !(updated.H ≈ model.H)
updated.coupling_parameters
```

Physical input changes reconstruct derived state; direct updates to derived `H`,
operators, or dressed states are rejected. Component updates must preserve their
names and order. Existing named gates are retained, so check their endpoint
constraints again after changing the idle model.
