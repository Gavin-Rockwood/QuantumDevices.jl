---
name: quantumdevices-models
description: Build QuantumDevices coupled devices, standalone recipes, and numerical models with explicit parameter overrides and state reductions.
---

Inspect existing components, parameter values, tensor order and truncations first. Construct only what is missing. Definitions are cheap; `build` generates operators, diagonalizes and may track interaction ramps. Preserve notebook setup/settings and do not build merely to inspect definitions.

## Complete small coupled model

Executable from the repository project. Frequencies/couplings below are illustrative E/h in GHz (time ns), not calibrated device values.

```julia
using QuantumDevices

device = QuantumDevice(:coupled)
q1 = add_component!(device, make_qubit(:q1; frequency=5.0))
q2 = add_component!(device, make_qubit(:q2; frequency=5.4))
add_interaction!(device, Interaction(:xx,
    0.01 * q1.operators["x"] * q2.operators["x"]))
recipe = Recipe(device, [q1, q2];
    states=Dict(q1 => 0:1, q2 => 0:1), interactions=:auto)
ν1 = q1.variables["ν"]
model = build(device, recipe; parameters=Dict(ν1 => 5.1))
H = make_hamiltonian(model)
ψ = state(model, (1, 0))
E = stateenergy(model, (1, 0))
@assert size(H) == (4, 4)
@assert model.parameters[ν1] == 5.1
@assert model.eigensystem.components == (:q1, :q2)
@assert isapprox(H * ψ, E * ψ; atol=1e-9)
```

`recipe` is an unregistered standalone definition; registration is optional (`add_recipe!(device, :pair, recipe)`). `model` is a NumericalModel with immutable resolved `parameters`, energy-ordered `spectrum`, labeled `states`/`energies`, numerical operators and cached complete products. `H` is a static QuantumToolbox QuantumObject in model construction coordinates; compilation does not rotate into the energy basis. `get_hamiltonian` instead returns a symbolic expression.

Components expose string-keyed `variables` (symbolic parameters), `parameters` (defaults), `operators` (symbolic declarations), and `operator_generators`. `interactions=:auto` captures applicable registered cross-input terms when constructing the recipe; default interactions are empty. Explicit `interactions=[Interaction(...)]` also works. Registration after construction does not alter that recipe. Scoped local additions use `recipe.inputs.q1.operators.x`; grouped inputs add the physical component segment.

## Choose the right change

- A new reference operating point: `build(device, recipe; parameters=Dict(variable => value))`. Do not mutate numerical caches to change physics. For saved definitions edit `draft(device.recipes[:pair])` and register/replace deliberately.
- Pure projected reduction: `reduced = Recipe(recipe; states=[(0,0), (1,0)])`, then `build(device, reduced; parameters=Dict(ν1 => 5.1))`. This builds the upstream model before selecting its dressed states; it does not avoid that upstream cost or repeat its couplings.
- Runtime coefficients in the existing fixed basis: read [make_hamiltonian](../quantumdevices-make-hamiltonian/SKILL.md). Runtime changes do not rebuild projections or relabel eigenstates.
- Two transmons plus a flux-tunable coupler: read [physical-device example](references/physical-device.md) only for that setup.
- Evolution: continue with [Controls](../quantumdevices-controls/SKILL.md). Persistence: [Device I/O](../quantumdevices-device-io/SKILL.md).

Selections and occupation labels are zero-based; ordinary Julia array indexing is one-based. Inspect `model.eigensystem.components` for actual tuple order (explicit grouped inputs precede uncovered components). `stateindex` returns a zero-based energy-order index, not a coordinate ket. Use `state(model, key)` instead of assuming `basis` matches eigenstates, including for nonzero-frequency Pauli qubits.

Projections capture a fixed construction basis. Complete operator words are multiplied before projection: `P' * A * B * P` need not equal `(P' * A * P) * (P' * B * P)`. Rebuild for a changed reference basis, generator dimensions or selections. Default H uses E/h; direct angular-frequency solvers need `scalar=2pi`, whereas built-in `evolve` supplies 2π once.

Read [quickstart](../../docs/src/getting_started/quickstart.md) for grouped composition, [models](../../docs/src/user_guide/models.md) for contracts, [parameters](../../docs/src/user_guide/parameters.md) for precedence, or [projection](../../docs/src/explanations/projection.md) when needed. Check cutoff/ramp convergence for a quantitative physical prediction; toy cutoffs are not convergence evidence.
