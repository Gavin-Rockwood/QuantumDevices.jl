# Physical and numerical conventions

## Hamiltonian units

The evolution convention is

```math
\frac{d}{dt}|\psi\rangle=-iH(t)|\psi\rangle,\qquad \hbar=1.
```

Hamiltonian coefficients therefore have units inverse to time. For time in ns,
an energy coefficient quoted as a frequency in GHz becomes `2π * frequency` in
rad/ns. The package does not apply that conversion implicitly.

`make_qubit(name, ν)` uses `H=ν Z/2`; its splitting is `|ν|` in angular-frequency
units. The first Pauli basis vector has Z eigenvalue +1. For positive `ν`, it is
not the lower-energy state. Energy-level product labels are determined by local
energy ordering and must not be substituted for computational basis indices.

`sine_pulse` uses `sin(2π*f*t + phase)`, with `f` in cycles/time and phase in radians.
`sine_squared_pulse` squares this sine, so its waveform repeats at twice `f`.
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

`unitary_infidelity` implements the first metric and clamps it to `[0,1]`.
It checks target unitarity, not actual unitarity. Validate solver accuracy and
unitarity when interpreting very small errors. The formula is phase insensitive.
It compares the complete retained space, not one state transfer and not a projected
computational subspace. It is not a general fidelity for dissipative channels.

## Numerical limits

Check local charge-cutoff convergence, retained-dimension convergence, state-tracking
step density, and evolution tolerances independently. A good result on one check
does not establish the others. Tracking confidence reflects continuity of sampled
matches; it is not a physical uncertainty estimate.
