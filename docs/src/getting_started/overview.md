# Installation and core concepts

QuantumDevices is currently developed from a repository checkout. In a directory
where you want to use it, open Julia and add the repository:

```julia
using Pkg
Pkg.add(url="https://github.com/Gavin-Rockwood/QuantumDevices.jl")
```

For development, activate the checkout instead:

```sh
julia --project=.
```

```julia
using Pkg
Pkg.instantiate()
using QuantumDevices
```

The package declares Julia 1.10 or newer; dependency resolution determines which
compatible versions can be installed. The CI build uses the current Julia release.
Use one project environment consistently when running scripts and notebooks.

## The workflow

1. Construct a **Component** with local operators, physical parameters, and a symbolic Hamiltonian.
2. Build a **DeviceModel** with components, their interaction, and retained dimensions.
3. Construct a **DeviceGate** with control parameters, an additional Hamiltonian, and duration.
4. Evaluate the gate Hamiltonian with `numerical(model, gate)` and use QuantumToolbox for evolution.
5. Calibrate selected scalar parameters or save the model and its named gates.

All public names live in `QuantumDevices`. Import QuantumToolbox explicitly for
operators, initial states, and solvers. No Circuits or Dynamics submodule is needed.

## Units and labels

Hamiltonians follow ℏ = 1. An input frequency in Hz must be converted to angular
frequency using `2π`; time and Hamiltonian coefficients must use reciprocal units.
Pulse constructors named `sine_pulse` and `sine_squared_pulse` instead take their
sine frequency in **cycles per time**. See [conventions](../explanations/conventions.md).

Tensor order follows component order. Dressed-state labels are tuples of
zero-based local **energy levels**, while Julia array indices are one-based.
These labels need not match a qubit's computational basis indices.

Continue with the [quickstart](quickstart.md).
