# Transmon resonator control

Build a coupled transmon–resonator model, attach shaped charge drives, and
simulate transfer from the dressed `|f,0⟩` state to `|g,1⟩`. This tutorial uses a
saved control set so you can focus on assembling and running the model.

## 1. Set up and load the controls

From the repository root, install the demo dependencies once:

```sh
julia --project=demo -e 'using Pkg; Pkg.instantiate()'
```

You can run the whole [Julia tutorial](https://github.com/Gavin-Rockwood/QuantumDevices.jl/blob/main/demo/mode3_controls.jl)
with `julia --project=demo demo/mode3_controls.jl`. It prints the final population
and saves two figures in a temporary directory; pass a directory as the last
argument to choose where to save them. To follow along interactively, start
`julia --project=demo` and use the steps below.

```@example mode3
using QuantumDevices, QuantumToolbox, CairoMakie, LinearAlgebra, SciMLBase
include(joinpath(pkgdir(QuantumDevices), "demo", "mode3_controls.jl"))
using .Mode3ControlsDemo
controls = load_controls()
config = controls["Main_Config"]
saved = controls["Stuff"]["op_drive_params"]["sb_f0g1"]
(; duration_ns=saved["pulse_time"], frequency_GHz=saved["freq_d"])
```

The [data file](https://github.com/Gavin-Rockwood/QuantumDevices.jl/blob/main/demo/data/mode3_controls.json)
contains model parameters, qubit pulses, and the `sb_f0g1` sideband. Time is in
ns; saved energies, frequencies, and amplitudes are in cycles/ns (GHz).
QuantumDevices Hamiltonians use angular units, so multiply these energies and
drive coefficients by `2π`.

## 2. Build the components and their interaction

Create a transmon in the saved charge basis and a resonator in its Fock basis.
Component names determine operator prefixes: `q` supplies `:q_charge`, and `r`
supplies `:r_a` and `:r_adag`.

```@example mode3
q = make_transmon("q", 2pi * config["E_C"], 2pi * config["E_J"],
    2 * Int(config["Nt_cut"]) + 1; ng=config["ng"])
r = make_resonator("r", 2pi * config["E_oscs"][1], Int(config["Nrs"][1]))
interaction = param(:g) * op(:q_charge) * (1im * (op(:r_a) - op(:r_adag)))
model = make_model([q, r], interaction, (; g=2pi * config["gs"][1]);
    truncation_dimensions=Dict(q => 4, r => 3), max_dimension=12)
size(model.H)
```

`make_model` includes each component Hamiltonian and adds the interaction. Here
we retain four transmon and three resonator levels for a quick construction
example. The transmon parent charge basis still uses the saved `Nt_cut=60`.
The interaction uses charge coupled to `i(a-a†)`, matching the saved model.

## 3. Attach a shaped drive

First reconstruct the envelope and carrier from one saved pulse record:

```@example mode3
envelope = saved_envelope(saved)
drive = drive_pulse(saved)
duration = Float64(saved["pulse_time"])
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
h(t)=2\pi\epsilon\,e(t)\sin\!\left[2\pi(f_{\mathrm d}+\mathrm{shift})t\right].
```

`shift` is a frequency offset in GHz. The helper combines a unit-height envelope
with this carrier using `GenericPulseFunction`; its callable receives the gate
duration as `p.duration`. The saved Gaussian records use `gaussian_pulse`;
their original `Guassian` spelling is retained in the JSON. The sideband uses
`ramped_flattop_pulse` with a custom bump source and the saved ramp time.

Inspect all the saved controls before evolving the system:

```@example mode3
plot_controls(controls)
```

The left panels show amplitude times envelope. The right panels show short
carrier segments, with Hamiltonian coefficients divided by `2π`.

## 4. Evolve the sideband transfer

The saved sideband was calibrated in a 10-by-10 retained space. Use
`make_mode3_model(controls)` to build that space and attach every saved control.
Its dimension keywords let you choose smaller spaces for exploratory runs,
but reproducing the saved transfer requires the full space.

```@example mode3
full_model = make_mode3_model(controls)
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
population transfer; the saved `accuracy` fields are metadata from the original
calibration. This run does not apply the save's damping parameters, other
sidebands, or chirp fit. The original controls came from a
[SuperconductingCavitiesDemo](https://github.com/Gavin-Rockwood/SuperconductingCircuits.jl)
simulation.
