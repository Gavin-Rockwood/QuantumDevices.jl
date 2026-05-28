# Overview

QuantumDevices.jl is a Julia package for building circuit-level models of quantum devices and simulating their driven dynamics. The current implementation is organized around two main workflows:

- Build circuit components, combine them into a composite Hamiltonian, and add named operators for later analysis.
- Use dynamics tools for propagators, Floquet bases and sweeps, drive construction, resonance searches, and drive calibration.

The package reexports the public APIs from `QuantumDevices.Circuits`, `QuantumDevices.Dynamics`, and `QuantumDevices.Utils`, so most examples start with:

```julia
using QuantumDevices
```

Some lower-level examples also use QuantumToolbox directly:

```julia
using QuantumToolbox
```

## Local Development Setup

From a checkout of the repository, use the root project for package development and tests:

```sh
julia --project=.
```

Inside Julia, instantiate the environment if needed:

```julia
using Pkg
Pkg.instantiate()
```

The documentation has its own project under `docs/`:

```sh
julia --project=docs docs/make.jl
npm --prefix docs run docs:build
```

## Core Concepts

Circuit components such as `Qubit`, `Resonator`, `Transmon`, `Fluxonium`, `FluxTunableTransmon`, and `SNAIL` store local Hamiltonians, eigenstates, dimensions, parameters, and loss operators. Constructors such as `init_transmon` and `init_resonator` create these components directly, while `init_components` provides a string-keyed registry for configuration-driven construction.

Circuits are built with `init_circuit(components, interactions; ...)`. Interactions and named operators use the typed operator-expression DSL:

```julia
op(:transmon, :n) * op(:resonator, :a)
coupling(0.01, op(:transmon, :n) * op(:resonator, :a); hc=true)
```

The symbols in `op(component, operator)` must match each component's `:name` parameter. Operators are resolved with `get_operator` and automatically tensored with identities on components that are not referenced.

Dynamics tools then operate on QuantumToolbox Hamiltonians, circuit Hamiltonians, and named circuit operators. The main pieces are:

- `propagator` and `Propagator` for time evolution operators.
- `get_floquet_basis` and `floquet_sweep` for periodic Hamiltonians.
- Envelope functions and drive parameter structs for building time-dependent drives.
- `find_resonance`, `calibrate_drive`, and `get_FLZ_flattop` for resonance and pulse-calibration workflows.

Start with the circuit guide if you are building Hamiltonians, or the dynamics guides if you already have Hamiltonians and drive operators.
