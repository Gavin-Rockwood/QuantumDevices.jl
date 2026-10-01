# Transmon resonator control

Build a coupled transmon–resonator model, attach shaped charge drives, and
simulate transfer from the dressed `|f,0⟩` state to `|g,1⟩`. This tutorial uses a
fixed parameter set defined directly in the Julia demo.

## 1. Set up and choose parameters

From the repository root, install the demo dependencies once:

```sh
julia --project=demo -e 'using Pkg; Pkg.instantiate()'
```

You can run the whole [Julia tutorial](https://github.com/Gavin-Rockwood/QuantumDevices.jl/blob/main/demo/mode3_controls.jl)
with `julia --project=demo demo/mode3_controls.jl`. It prints the final population
and saves a model bundle and two figures in a temporary directory; pass an
output directory as the last argument to choose where to save them. To follow
along interactively, start `julia --project=demo` and use the steps below.

```@example mode3
using QuantumDevices, QuantumToolbox, CairoMakie, LinearAlgebra, SciMLBase
include(joinpath(pkgdir(QuantumDevices), "demo", "mode3_controls.jl"))
using .Mode3ControlsDemo
params = TARGET
settings = SIDEBAND_SETTINGS
(; model_parameters=params, sideband_parameters=settings)
```

`TARGET` defines the transmon energies, resonator frequency, coupling, and basis
sizes. `QUBIT_SETTINGS` defines two Gaussian controls, and `SIDEBAND_SETTINGS`
defines the sideband. All values are declared in `demo/mode3_controls.jl`.
Time is in ns; energies, frequencies, and amplitudes are in cycles/ns (GHz).
QuantumDevices Hamiltonians use angular units, so multiply energies and drive
coefficients by `2π`.

## 2. Build the components and their interaction

Create a transmon in the charge basis and a resonator in its Fock basis.
Component names determine operator prefixes: `q` supplies `:q_charge`, and `r`
supplies `:r_a` and `:r_adag`.

```@example mode3
q = make_transmon("q", 2pi * params.EC, 2pi * params.EJ,
    2params.n_cutoff + 1; ng=params.ng)
r = make_resonator("r", 2pi * params.resonator_frequency, params.resonator_levels)
interaction = param(:g) * op(:q_charge) * (1im * (op(:r_a) - op(:r_adag)))
model = make_model([q, r], interaction, (; g=2pi * params.g);
    truncation_dimensions=Dict(q => 4, r => 3), max_dimension=12)
size(model.H)
```

`make_model` includes each component Hamiltonian and adds the interaction. Here
we retain four transmon and three resonator levels for a quick construction
example. The transmon parent charge basis uses `n_cutoff=60`.
The interaction couples transmon charge to the resonator quadrature `i(a-a†)`.

## 3. Attach a shaped drive

Build the sideband envelope and its carrier from the declared pulse parameters:

```@example mode3
envelope = sideband_envelope(settings)
drive = drive_pulse(settings, envelope)
duration = settings.duration
gate = DeviceGate((; drive), param(:drive) * op(:q_charge), duration)
model.gates[:sb_f0g1] = gate
(; envelope_at_center=pulse_value(envelope, duration / 2, duration),
   duration_ns=gate.duration)
```

A `DeviceGate` pairs parameter values with a symbolic drive Hamiltonian and a
duration. `param(:drive)` reads the pulse stored under `drive`.
`pulse_value(pulse, t, duration)` evaluates a pulse at a physical time.

The carrier coefficient is

```math
h(t)=2\pi\epsilon\,e(t)\sin(2\pi\nu t).
```

`settings.frequency` is the carrier frequency `ν` in GHz. The helper combines
a unit-height envelope with this carrier using `GenericPulseFunction`; its
callable receives the gate duration as `p.duration`. The qubit controls use
`gaussian_pulse`, and the sideband uses `ramped_flattop_pulse` with a custom bump
source and the declared ramp time.

Build the complete model and inspect its three controls before evolving it:

```@example mode3
full_model = make_mode3_model()
plot_controls(full_model)
```

The left panels show amplitude times envelope. The right panels show short
carrier segments, with Hamiltonian coefficients divided by `2π`.

## 4. Evolve the sideband transfer

Use ten retained levels per component to resolve the sideband transfer with
these pulse parameters. `make_mode3_model()` builds that space and attaches the
two Gaussian drives and the sideband. Its dimension keywords let you choose
smaller spaces for exploratory runs.

```@example mode3
gate = full_model.gates[:sb_f0g1]
initial = full_model.states[(2, 0)]
target = full_model.states[(0, 1)]
times = range(0, gate.duration; length=201)
result = sesolve(numerical(full_model, gate), initial, times;
    progress_bar=false, abstol=1e-8, reltol=1e-8)
@assert SciMLBase.successful_retcode(result.retcode)
population = [abs2(dot(target, state)) for state in result.states]
@assert last(population) > 0.99
(; final_population=last(population), dressed_gap_GHz=real(sideband_gap(full_model)))
```

State labels are zero-based and follow the component order `[q, r]`:
`(2, 0)` is `|f,0⟩`, and `(0, 1)` is `|g,1⟩`. `model.states` supplies dressed
states. `numerical(model, gate)` combines the model and control into the
Hamiltonian consumed by QuantumToolbox's `sesolve`. The squared overlap with
the target gives its population at each sample time. The script's
`sideband_population` helper wraps these same steps.

```@example mode3
transfer = Figure(size=(700, 380))
ax = Axis(transfer[1, 1]; xlabel="Time (ns)", ylabel="|g,1⟩ population",
    title="f0g1 sideband transfer")
lines!(ax, times, population; color=:purple, linewidth=3)
ylims!(ax, 0, 1.05)
transfer
```

You should obtain a final target population above 0.99. This is closed-system
population transfer with the declared pulse parameters; no dissipation or
calibration step is applied.

## 5. Save and reload the model

Save the instantiated model, including all three named gates, as a model bundle:

```@example mode3
bundle_path = QuantumDevices.save(joinpath(mktempdir(), "transmon_resonator"), full_model)
restored = QuantumDevices.load(bundle_path)
@assert restored.H ≈ full_model.H
@assert Set(keys(restored.gates)) == Set(keys(full_model.gates))
@assert pulse_value(restored.gates[:sb_f0g1].parameters.drive, duration / 2, duration) ≈
    pulse_value(full_model.gates[:sb_f0g1].parameters.drive, duration / 2, duration)
println(join(sort(readdir(bundle_path)), ", "))
```

`QuantumDevices.save` creates a new directory and refuses an existing destination.
The standalone demo saves this bundle in `<output_directory>/model/`, alongside
`controls.png` and `transfer.png`, and prints the path.

In a fresh Julia session, include the demo module before loading the bundle so
its custom pulse callables are defined:

```julia
using QuantumDevices
include(joinpath(pkgdir(QuantumDevices), "demo", "mode3_controls.jl"))
restored = QuantumDevices.load("output_directory/model")
```

Including the file defines its module without running the simulation or writing
files. See [Saving and loading model bundles](../user_guide/persistence.md) for
the bundle layout.
