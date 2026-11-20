# Physical and numerical conventions

## Hamiltonian units

All Hamiltonians and dimensional energy parameters use frequency units, `E/h`,
in **cycles per unit time**, not angular frequency. With time in ns, use GHz.
The evolution convention is

```math
\frac{d}{dt}|\psi\rangle=-2\pi iH(t)|\psi\rangle.
```

Pass frequencies directly to component constructors and gate coefficients.
`make_hamiltonian` preserves frequency units by default, including the model's dressed eigenvalues.
Built-in `evolve` supplies the `2π` conversion internally. When calling
QuantumToolbox `sesolve` or `mesolve` directly, use
`make_hamiltonian(model, gate; scalar=2pi)` to fold `2π` into the matrices, or supply
`2π * H` yourself;
Do not add another `2π` factor when using `evolve`.

`make_qubit(name, ν)` uses `H=ν Z/2`; its splitting is `|ν|` in cycles per time.
The first Pauli basis vector has Z eigenvalue +1. For positive `ν`, it is
not the lower-energy state. Energy-level product labels are determined by local
energy ordering and must not be substituted for computational basis indices.

`Pulse(shape, f; phase=0)` mixes an I+iQ envelope as `I*cos(2π*f*t+phase)-Q*sin(2π*f*t+phase)`, with f in cycles/time and phase in radians.
Its clock defaults to pulse-local time. `SineSquared(; duration=T)` is a single envelope lobe over the pulse duration.
For flux-tunable transmons, `phi` is flux divided by the flux quantum and `ng` is
charge offset in the constructor's Cooper-pair convention.

## Tensor and frame conventions

Components `[a, b]` are ordered `a ⊗ b`. Gate targets and initial states must use
that order and the model's retained dimensions. There is no automatic rotating
frame or rotating-wave approximation. Explicit resonance analysis is available through
`find_resonance`; it does not change the gate evolution frame. A drive is a Hamiltonian
coefficient multiplying its specified operator; its amplitude need not equal the
observable Rabi angular frequency without an operator normalization convention.

## Gate metrics

For target and realized unitaries on a d-dimensional retained space,

```math
1-F_{\mathrm{pro}}=1-\frac{|\operatorname{Tr}(U_{\mathrm{target}}^\dagger U)|^2}{d^2},
\qquad
1-F_{\mathrm{avg}}=\frac{d}{d+1}(1-F_{\mathrm{pro}}).
```

The current `gate_infidelity(T,M)` uses the leakage-sensitive average loss

```math
1-\frac{\operatorname{Tr}(M^\dagger M)+|\operatorname{Tr}(T^\dagger M)|^2}{d(d+1)}.
```

The target is unitary and input/output order is explicit. Projected columns are
never normalized: their missing norm measures leakage. Relative phases count;
a common global phase cancels. `gate_matrix(evolution)` returns the physical
projected propagator, and `idle_phases(evolution)` reports reference idle phases
without applying a frame correction. No automatic interaction-frame or virtual-Z
compensation is performed. These are closed-system scores, not dissipative-channel fidelities.

## Numerical limits

Check local charge-cutoff convergence, retained-dimension convergence, state-tracking
step density, and evolution tolerances independently. A good result on one check
does not establish the others. Tracking confidence reflects continuity of sampled
matches; it is not a physical uncertainty estimate.
