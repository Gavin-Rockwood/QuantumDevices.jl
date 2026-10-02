# Components

A component owns a parent Hilbert space, operators in a consistent local basis,
physical parameters, and a symbolic Hamiltonian. The built-in constructors are
`make_qubit`, `make_resonator`, `make_transmon`, and `make_tunable_transmon`.

```@example components
using QuantumDevices, QuantumToolbox
q = make_qubit("q", 1.0)
r = make_resonator("r", 2.5, 4)
t = make_transmon("t", 0.2, 5.0, 9)
ft = make_tunable_transmon("ft", 0.2, 5.0, 4.0, 9; phi=0.2)
@assert keys(t.operators) == (:jump, :charge)
@assert r.operators.n ≈ r.operators.adag * r.operators.a
(q.dimension, r.dimension, t.dimension, ft.dimension)
```

| Component | Parameters | Local operators | Parent basis |
|:--|:--|:--|:--|
| Qubit | `ν` | `x`, `y`, `z`, `p`, `m` | Two-level Pauli basis |
| Resonator | `frequency` | `a`, `adag`, `n` | Fock basis, dimension is the number of included states |
| Transmon | `EC`, `EJ`, `ng` | `jump`, `charge` | Energy basis obtained from a finite charge basis |
| Tunable transmon | `EC`, `EJ1`, `EJ2`, `ng`, `phi`; derived `EJmax`, `d` | `jump`, `charge` | Initial energy basis at the construction flux |

Use `numerical(component, expression)` with the component's local operator names:

```@example components
H = numerical(t, t.hamiltonian)
charge = numerical(t, op(:charge))
@assert charge ≈ t.operators.charge
size(H)
```

## Parent dimensions and convergence

Resonator dimensions must be positive integers. Increase the Fock cutoff until
the occupations and dynamics of interest converge. Transmon dimensions must be
positive odd integers: `dimension = 2ncut + 1`.
The charge basis covers `-ncut:ncut`. This cutoff controls convergence of the
local spectrum; the model's retained dimension is a separate approximation.
Increase the parent cutoff first, then check retained-model predictions.

The constructor transforms charge and Josephson operators into its initial energy
basis. Changing physical construction parameters with `setpath` rebuilds that
basis. A time-dependent parameter in a gate acts in the already constructed basis;
it does not diagonalize a new basis at each time sample.

## Custom components

Use `Component(name, parameters, operators, hamiltonian, type, dimension)` for
additional device types. Operators must match the parent dimension, and local names
must match the symbolic Hamiltonian. See [extensions](../development/extensions.md)
and the [component API](../reference/models.md).
