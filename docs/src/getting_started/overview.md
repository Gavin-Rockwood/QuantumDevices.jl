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

The package requires Julia 1.11 or newer in the Julia 1.x series. CI tests Julia
1.11 and the current Julia release. Dependency resolution selects compatible
package versions.
Use one project environment consistently when running scripts and notebooks.

## The workflow

1. Construct components and register them with a `QuantumDevice`.
2. Define a `Recipe` and build its numerical model.
3. Define a `GateRecipe` with controls and optional interaction terms.
4. Build the gate or compile it against a compatible cached model, then evolve it.
5. Calibrate selected control parameters when needed.

All public names live in `QuantumDevices`. Import QuantumToolbox explicitly for
operators, initial states, and solvers. No Circuits or Dynamics submodule is needed.

## Units and labels

Hamiltonians and energy parameters use frequency units, `E/h`, in **cycles per
time** (GHz for time in ns). Pass frequency values directly to constructors.
`make_hamiltonian` preserves them by default; built-in `evolve` supplies the `2π`
conversion. For direct QuantumToolbox evolution, use `make_hamiltonian(...; scalar=2pi)`.
Carrier frequencies use cycles per time and phases use radians. See [conventions](../explanations/conventions.md).

Tensor order follows component order. Dressed-state labels are tuples of
zero-based local **energy levels**, while Julia array indices are one-based.
These labels need not match a qubit's computational basis indices.

Continue with the [quickstart](quickstart.md).
