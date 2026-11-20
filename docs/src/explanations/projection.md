# Projection and model bases

Let P map a retained basis into a component's parent basis. A retained operator is

```math
A_{\mathrm{ret}}=P^\dagger A P.
```

For a product, the physically compiled expression is

```math
(AB)_{\mathrm{ret}}=P^\dagger ABP,
```

which generally differs from `(P† A P)(P† B P)`: the latter inserts `PP†` between
factors and removes parent-space excursions. QuantumDevices groups each ordered
operator word by component, multiplies in its parent basis, projects once, and
then tensors the retained factors.

```@example projection
using QuantumDevices, QuantumToolbox, LinearAlgebra
A = QuantumObject([0.0 1 0; 1 0 1; 0 1 0])
operators, operator_generators = component_operators(:c;
    x=OperatorGenerator((dimension;kwargs...)->A,3))
c = Component(:c; variables=Dict(), parameters=Dict(), operators,
    operator_generators, hamiltonian=val(0)*operators["x"])
device = QuantumDevice(:projection)
add_component!(device,c)
recipe = Recipe(device,[c];states=Dict(c=>0:1),dressing=false)
gate = build(device,GateRecipe(recipe;duration=1.,
    terms=[Interaction(:square,recipe.inputs.c.operators.x^2)]))
compiled = make_hamiltonian(gate)
naive = gate.model.operators[OperatorRef(:c,:x)]^2
@assert !(compiled ≈ naive)
(Matrix(compiled.data), Matrix(naive.data))
```

This distinction applies to both idle and gate expressions. Individual
`model.operators` are useful for measurements and inspection, but multiplying
them is not a substitute for compiling a parent-space symbolic operator word.

Built-in transmon leaf operators stay in their charge construction basis. A
composite explicitly selects eigenvector columns from each upstream model;
these columns form `P`, and their tensor product defines its new construction
basis. The composite Hamiltonian and complete projected products stay in that
basis after dressing/diagonalization. Its exposed eigenstates use the same
coordinates, so `H * state = energy * state` holds directly.

A pure reduction from one upstream coupled recipe uses selected dressed
eigenvectors as its construction basis. It is naturally diagonal, preserving
inherited labels and diagnostics without another solve. Changing physical inputs
may change the selected upstream eigenvectors; compare physical observables or
use explicit transports when comparing independently rebuilt bases.
