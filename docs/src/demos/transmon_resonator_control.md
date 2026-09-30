# Transmon resonator control

This demo rebuilds the qubit Gaussian pulses and the `sb_f0g1` sideband from a
saved Mode3 control set. The [control data](https://github.com/Gavin-Rockwood/QuantumDevices.jl/blob/main/demo/data/mode3_controls.json)
and [Julia implementation](https://github.com/Gavin-Rockwood/QuantumDevices.jl/blob/main/demo/mode3_controls.jl)
are in `demo/` so the same code can be run outside the documentation site. The
data file contains the qubit records, `sb_f0g1`, and the model parameters from
the supplied save; it does not include the other sidebands or chirp fit.

The save came from a [SuperconductingCavitiesDemo](https://github.com/Gavin-Rockwood/SuperconductingCircuits.jl)
simulation. Its `Guassian` spelling is preserved in the data. The envelope
shapes are a Gaussian centered at `mu` and a bump-ramp with the saved `k` and
`ramp_time`. The latter uses this package's `ramped_flattop_pulse` with a bump
source. `shift` is a **frequency offset**, added to `freq_d`; it is not an
envelope offset. Time is in ns and saved frequencies and amplitudes are in
cycles/ns (GHz). The Hamiltonian coefficient is
`2π epsilon × envelope(t) × sin(2π (freq_d + shift) t)` in rad/ns.

```@example mode3
using QuantumDevices, CairoMakie
include(joinpath(pkgdir(QuantumDevices), "demo", "mode3_controls.jl"))
using .Mode3ControlsDemo
controls = load_controls()
fig = plot_controls(controls)
fig
```

The top row shows the `q_ge_0` and `q_ef_0` pulses. The lower row shows the
`sb_f0g1` envelope and a short section of its carrier. The carrier panels show
the Hamiltonian coefficient divided by `2π`.

```@example mode3
model = make_mode3_model(controls; transmon_levels=4, resonator_levels=3)
@assert length(model.gates) == length(controls["Stuff"]["op_drive_params"])
@assert abs(real(sideband_gap(model)) - abs(controls["Stuff"]["op_drive_params"]["sb_f0g1"]["freq_d"])) < 1e-3
(; qubit_gate=model.gates[:q_ge_0].duration,
   sideband_gate=model.gates[:sb_f0g1].duration,
   dressed_gap_GHz=real(sideband_gap(model)))
```

The model uses the saved transmon and resonator energies and their charge–mode
coupling. It drives the transmon charge operator for both qubit and sideband
controls. The example retains four transmon and three resonator levels to keep
the documentation build quick. Omit both dimension keywords to use the saved
10-by-10 retained space. Its charge-basis cutoff remains `Nt_cut=60` in either
case. The saved `accuracy` fields report results from the original calibration;
they are metadata, not fidelities measured by this example.

## Saved f0g1 transfer

The full retained space is needed to replay the saved sideband calibration.
Starting in the dressed `|f,0⟩` state, apply the saved `sb_f0g1` gate and track
the dressed `|g,1⟩` population:

```@example mode3
full_model = make_mode3_model(controls)
times, population = sideband_population(full_model; samples=201, abstol=1e-8, reltol=1e-8)
@assert last(population) > 0.99
transfer = Figure(size=(700, 380))
ax = Axis(transfer[1, 1]; xlabel="Time (ns)", ylabel="|g,1⟩ population",
    title="Saved f0g1 sideband")
lines!(ax, times, population; color=:purple, linewidth=3)
ylims!(ax, 0, 1.05)
transfer
```

This is coherent evolution under the saved Hamiltonian and drive. The source
save also lists damping parameters, but this replay does not apply dissipation
or the other saved pulse sequences. The numerical population is separate from
the save's `accuracy` field.
