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
gate = DeviceGate((drive=Pulse(RampedFlattop(0.2); amplitude=0.1, duration=1.0),),
                  param(:drive) * op(:q_x), 1.0)
H = numerical(model, gate)
H(0.5)
```

Hamiltonians and energy parameters use frequency units (`E/h`, cycles per unit
time): GHz for time in ns. `numerical` preserves these units. `get_unitary`
applies `2π` at the solver call; direct QuantumToolbox `sesolve` or `mesolve`
calls must use `2π * H`.
The current API replaces the old Circuits/Dynamics modules; start with the
[package overview](https://gavin-rockwood.github.io/QuantumDevices.jl/dev/getting_started/overview)
when upgrading existing code.

Run tests with `julia --project=. -e 'using Pkg; Pkg.test()'`. See
[development instructions](docs/src/development/contributing.md) for building the
site and extending components or controls.
