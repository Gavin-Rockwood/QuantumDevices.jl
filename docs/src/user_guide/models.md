# Coupled models and truncation

A Recipe selects registered components, their parent dimensions and the zero-based
energy states retained at each input. `build(device, recipe)` produces a NumericalModel.

```@example models
using QuantumDevices
qdevice = QuantumDevice(:coupled)
q = add_component!(qdevice, make_qubit(:q; frequency=1.0))
r = add_component!(qdevice, make_resonator(:r; frequency=2.0, dimension=3))
recipe = Recipe(qdevice, [q,r]; states=Dict(q=>0:1,r=>0:1), dressing=false)
add_interaction!(recipe, Interaction(:coupling,
    0.01 * recipe.inputs.q.operators.x *
    (recipe.inputs.r.operators.a + recipe.inputs.r.operators.a')))
model = build(qdevice, recipe)
@assert size(model.H) == (4,4)
size(model.H)
```

Selections are energy indices, while labels identify physical states. Set `initial_dims`
for automatic component input generators; retained states must fit their parent dimensions.
`max_dimension` bounds the constructed tensor dimension, not every subsequent solver allocation.

## Construction coordinates

For recipe-based models, `model.H` and `make_hamiltonian(model)` use the same construction
coordinates as `model.operators`, complete cached products, and `model.states`.
`model.H` is stored once at build, including the default numerical Hermiticity
cleanup. Default `make_hamiltonian(model)` returns an independent copy of that reference;
custom scaling or cleanup compiles raw cached products. Mutating operator/cache
arrays does not update stored H or its eigensystem; rebuild the model for a new
reference. No duplicate diagonalizing unitary is stored.
`state(model, key)` returns that same physical eigenvector, while `stateindex`
is only its zero-based energy-order index. Diagonalization provides energies and
labels without rotating model operators. Explicit downstream selections project
onto upstream eigenvectors and define a new construction basis; a pure reduction
is naturally diagonal. Builds still eagerly supply spectra and labels.

## Model properties

`model.parameters` is an immutable resolved parameter snapshot; `model.operators`
holds projected numerical operators, and `model.hamiltonian` is symbolic.
`model.H` is the numerical reference, `model.states` the labeled eigenvectors and
`model.spectrum` their energies in frequency units.
Use `state(model,label)`, `stateenergy(model,label)` and `stateindex(model,label)`
for physical label access. Inspect label diagnostics before trusting continuity near crossings.

```@example models
@assert length(model.spectrum) == length(model.states)
model.spectrum
```

A complete operator word is multiplied in the parent basis **before** projection.
Consequently, `make_hamiltonian` of a symbolic charge operator squared can differ from squaring
the corresponding projected numerical charge operator. Use symbolic gate/model expressions when the distinction
matters. See [projection](../explanations/projection.md).

Use `make_hamiltonian(model)` to compile the model Hamiltonian, or
`make_hamiltonian(model, [parameter])` for runtime free parameters. Both forms accept
`scalar=2pi` for angular-frequency units. See the [quickstart](../getting_started/quickstart.md).

## Changing inputs

Change physical inputs through a new build or an editable draft of a captured recipe:

```@example models
updated = build(qdevice, recipe; parameters=Dict(q.variables["ν"]=>1.1))
@assert model.parameters[q.variables["ν"]] == 1.0
@assert updated.parameters[q.variables["ν"]] == 1.1
updated.spectrum
```

Do not edit derived numerical Hamiltonians or cached operators to represent new physics.
A rebuild refreshes the reference basis and labels. Gate definitions live on the device,
with their own captured model recipe; rebuild the appropriate gate context explicitly.

## Follow a small avoided crossing

A rebuild changes the physical operating point and recomputes its basis. For this
four-dimensional qubit–resonator example, vary the qubit frequency through the
resonator frequency and plot the two single-excitation energy branches. These
are energy-ordered curves; an energy index alone does not identify a continuous
physical state through the avoided crossing.

```@example models
using CairoMakie
frequencies = range(1.8,2.2;length=11)
spectra = [collect(build(qdevice,recipe;parameters=Dict(q.variables["ν"]=>ν)).spectrum)
    for ν in frequencies]
fig = Figure(size=(700,320))
ax = Axis(fig[1,1];xlabel="Qubit frequency",ylabel="Energy above ground")
for row in 2:3
    lines!(ax,frequencies,[E[row]-E[1] for E in spectra];label="Energy branch $(row-1)")
end
axislegend(ax; position=:lt)
fig
```

For physical state continuity, inspect the recipe's interaction-ramp labels and
[tracking diagnostics](tracking.md); increase cutoffs and sampling near small gaps.
For a full driven-control application, continue to the
[transmon–resonator tutorial](../tutorials/transmon_resonator_control.md).
