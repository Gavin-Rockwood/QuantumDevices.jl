# QuantumDevices.jl

[![Build Status](https://github.com/Gavin-Rockwood/QuantumDevices.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Gavin-Rockwood/QuantumDevices.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://github.com/Gavin-Rockwood/QuantumDevices.jl/actions/workflows/Documentation.yml/badge.svg?branch=dev)](https://github.com/Gavin-Rockwood/QuantumDevices.jl/actions/workflows/Documentation.yml)

Build symbolic quantum-device Hamiltonians, project coupled models into retained
component bases, track dressed states, and simulate or calibrate time-dependent gates.

- [Documentation](https://gavin-rockwood.github.io/QuantumDevices.jl/dev/)
- [Quickstart](https://gavin-rockwood.github.io/QuantumDevices.jl/dev/getting_started/quickstart)
- [API reference](https://gavin-rockwood.github.io/QuantumDevices.jl/dev/resources/api)

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Gavin-Rockwood/QuantumDevices.jl")
```

From a checkout, activate the root environment with `julia --project=.` and run
`Pkg.instantiate()`.

## A model and shaped control

```julia
using QuantumDevices
q = make_qubit("q", 1.0)
model = make_model([q], val(0), (;))
gate = DeviceGate((drive=ramped_flattop_pulse(0.1, 0.2),),
                  param(:drive) * op(:q_x), 1.0)
H = numerical(model, gate)
H(0.5)
```

Hamiltonian coefficients use energy/ℏ units reciprocal to simulation time.
The current API replaces the old Circuits/Dynamics modules; see the documentation's
migration guide for supported equivalents.

Run tests with `julia --project=. -e 'using Pkg; Pkg.test()'`. See
[development instructions](docs/src/development/contributing.md) for building the
site and extending components or controls.
