# Physical and numerical conventions

## Hamiltonian units

All Hamiltonians and dimensional energy parameters use frequency units, `E/h`,
in **cycles per unit time**, not angular frequency. With time in ns, use GHz.
The evolution convention is

```math
\frac{d}{dt}|\psi\rangle=-2\pi iH(t)|\psi\rangle.
```

Pass frequencies directly to component constructors and gate coefficients.
`numerical` preserves frequency units by default, including the model's dressed eigenvalues.
`get_unitary` explicitly supplies `2π * H` to QuantumToolbox. When calling
QuantumToolbox `sesolve` or `mesolve` directly, use
`numerical(model, gate; scalar=2pi)` to fold `2π` into the matrices, or supply
`2π * H` yourself;
do not also multiply constructor inputs or drive amplitudes by `2π`.

`make_qubit(name, ν)` uses `H=ν Z/2`; its splitting is `|ν|` in cycles per time.
The first Pauli basis vector has Z eigenvalue +1. For positive `ν`, it is
not the lower-energy state. Energy-level product labels are determined by local
energy ordering and must not be substituted for computational basis indices.

`SineCarrier` uses `sin(2π*f*t + phase)`, with `f` in cycles/time and phase in radians.
Its clock defaults to pulse-local time. `SineSquared()` is a single envelope lobe over the pulse duration.
For flux-tunable transmons, `phi` is flux divided by the flux quantum and `ng` is
charge offset in the constructor's Cooper-pair convention.

## Tensor and frame conventions

Components `[a, b]` are ordered `a ⊗ b`. Gate targets and initial states must use
that order and the model's retained dimensions. There is no automatic rotating
frame, resonance search, or rotating-wave approximation. A drive is a Hamiltonian
coefficient multiplying its specified operator; its amplitude need not equal the
observable Rabi angular frequency without an operator normalization convention.

## Gate metrics

For target and realized unitaries on a d-dimensional retained space,

```math
1-F_{\mathrm{pro}}=1-\frac{|\operatorname{Tr}(U_{\mathrm{target}}^\dagger U)|^2}{d^2},
\qquad
1-F_{\mathrm{avg}}=\frac{d}{d+1}(1-F_{\mathrm{pro}}).
```

`unitary_fidelity` returns `F_pro`; `unitary_infidelity` returns `1-F_pro`,
with results clamped to `[0,1]`.
It checks target unitarity, not actual unitarity. Validate solver accuracy and
unitarity when interpreting very small errors. The formula ignores global phase
but includes relative phases. For a projected gate matrix `M`, the same overlap
score retains leakage; the average-fidelity conversion above applies to unitaries.

With `include_phases=false`, the score is

```math
F_{\mathrm{prob}}=\frac{1}{d}\sum_j\left(\sum_i |T_{ij}|\,|M_{ij}|\right)^2.
```

This compares each column's transition probabilities without renormalizing
projected columns. For an X or permutation target it is the mean correct-transfer
probability. Mean leakage is `1 - sum(abs2, M)/d` for normalized evolved inputs
and an orthonormal output basis. These are closed-system gate scores, not general
fidelities for dissipative channels.

`get_gate_matrix` and `CalibrationProblem` support selected input/output bases.
Their default `frame=:lab` preserves lab phases; explicit `frame=:interaction`
computes `C' * exp(2π*im*H_idle*T) * Ψ(T)`. It changes the final comparison,
not the integration or its computational cost.

## Numerical limits

Check local charge-cutoff convergence, retained-dimension convergence, state-tracking
step density, and evolution tolerances independently. A good result on one check
does not establish the others. Tracking confidence reflects continuity of sampled
matches; it is not a physical uncertainty estimate.
