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

## Build a model

```julia
using QuantumDevices
using QuantumToolbox

device = QuantumDevice(:two_qubits)
q1 = add_component!(device, make_qubit(:q1; frequency=5.0))
q2 = add_component!(device, make_qubit(:q2; frequency=6.2))
recipe = Recipe(device, [q1, q2]; states=Dict(q1 => 0:1, q2 => 0:1))
add_interaction!(recipe, Interaction(:coupling,
    0.02 * recipe.inputs.q1.operators.x * recipe.inputs.q2.operators.x))
model = build(device, recipe)
H = make_hamiltonian(model) # static QuantumToolbox QuantumObject in the construction basis
```

Physical coefficients use Symbolics.jl directly: `q1.variables["ν"]` is the owned
symbolic frequency and `q1.parameters["ν"]` is its numerical idle value.
Use the variable as an override key, for example
`build(device, recipe; parameters=Dict(q1.variables["ν"]=>5.2))`.

`Recipe` defines the inputs, retained states, and local interactions without
numerical work. `build` diagonalizes the model; recipes can also consume other
recipes as named inputs. Device dictionaries support lookups such as
`device.components[:q1]` and `device.recipes[:target]` after registration. See the [quickstart](docs/src/getting_started/quickstart.md)
for a nested transmon example. `Leaf`, `Combine`, and `Input` remain lower-level
compatibility interfaces. Built models expose `model.states` and `model.energies` tracked arrays with
interaction-ramp occupation-tuple lookup, for example `model.energies(1, 0)`;
see [state labeling](docs/src/user_guide/state_labeling.md) for accessors and diagnostics.

Registration captures an immutable recipe tree. Create an editable trial with
`trial = draft(device.recipes[:target])`, change `trial.parameters`, and build it.
Use `replace_recipe!(device, :target, trial)` to capture an explicit replacement.
Saved recipes and built model provenance keep their captured definitions. Callables
may capture mutable external state, which remains the caller's responsibility.

Hamiltonians and energy parameters use frequency units (`E/h`, cycles per unit
time): GHz for time in ns. `numerical` preserves these units. Use
`make_hamiltonian(model, [parameter])` for a `QobjEvo` with runtime parameter functions
(`auto_time_dependence=false` reads scalar runtime values). Supply `2π * H`
to QuantumToolbox time-evolution solvers.

## Share a device

```julia
save_device("devices/two_qubits", device)
restored = load_device("devices/two_qubits")
```

Register working gates with `add_gate!` at named or nested variant paths; attach
saved evolutions explicitly with `attach_artifact!`. Complete folders reuse recipe
definitions and retain independent historical calibration snapshots. Loading does
no model build or simulation. Custom bundled source requires `allow_code=true`;
optional native optimizer checkpoints have a separate trust opt-in. See the
[persistence guide](docs/src/user_guide/persistence.md).

Run tests with `julia --project=. -e 'using Pkg; Pkg.test()'`. See
[development instructions](docs/src/development/contributing.md) for building the
site and extending components or controls.


The documentation includes worked [gate evolution](docs/src/user_guide/gates.md)
and full [transmon–resonator](docs/examples/tutorials/transmon_resonator_control.jl)
and [tunable-coupler](docs/examples/tutorials/tunable_coupler_control.jl) control
workflows. Running the full tutorial scripts performs substantial numerical work;
the site builds use only lightweight illustrative examples.
