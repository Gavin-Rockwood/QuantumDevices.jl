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
c = Component("c", (;), (; x=A), val(0) * op(:x), "custom", 3)
model = make_model([c], val(0), (;); truncation_dimensions=Dict(c=>2))
gate = DeviceGate((;), op(:c_x)^2, 1.0)
compiled = numerical(model, gate)
naive = model.operators.c_x^2
@assert !(compiled ≈ naive)
(Matrix(compiled.data), Matrix(naive.data))
```

This distinction applies to both idle and gate expressions. Individual
`model.operators` are useful for measurements and inspection, but multiplying
them is not a substitute for compiling a parent-space symbolic operator word.

Built-in transmon parent operators are constructed in the energy basis of the
initial physical parameters. Retaining its first d vectors selects its lowest
initial local energy levels. A physical reconstruction can change that basis;
compare physical observables and spectra rather than interpreting entries of
changed-basis matrices as fixed coordinates.
