# Gates and evolution

A `DeviceGate` contains parameters and an additional Hamiltonian. Its duration
is inferred from the latest pulse end unless supplied explicitly.
Parameters must be finite numbers or pulses; raw functions are rejected. Use
`Pulse(Envelope(f, parameters); duration)` for custom `f(p, t, duration)` envelopes.

```@example gates
using QuantumDevices, QuantumToolbox, SciMLBase, LinearAlgebra
q = make_qubit("q", 0.0)
model = make_model([q], val(0), (;))
gate = DeviceGate((; drive=Pulse(Constant(); amplitude=0.25, duration=1.0)),
                  param(:drive) * op(:q_x))
U = get_unitary(model, gate; abstol=1e-9, reltol=1e-9)
@assert unitary_infidelity(-im * sigmax(), U) < 1e-7
unitary_infidelity(-im * sigmax(), U)
```

The ideal zero-idle example isolates an X rotation. In a nonzero idle model, lab
frame detuning and drift remain present unless explicitly modeled or compensated.
There is no automatic rotating frame or rotating-wave approximation.

For state evolution, use the same Hamiltonian with QuantumToolbox:

```@example gates
result = sesolve(model, gate, basis(2, 0), range(0, 1; length=21); dense=true, progress_bar=false)
@assert SciMLBase.successful_retcode(result.retcode)
@assert abs2(dot(basis(2, 1), result.states[end])) > 0.999
length(result.states)
```

`numerical(...; scalar=2pi)` folds the frequency-to-angular-frequency conversion
into the operator matrices before building `QobjEvo`, avoiding lazy outer scaling.
Use `numerical(model, gate; scalar=2pi, dense=true)` to materialize every static
and time-dependent term's matrix as a dense `Matrix` before constructing the
Hamiltonian. Conversion happens once; coefficient callbacks only evaluate the
signals. The device-aware solver accepts the same option:
`sesolve(model, gate, initial, times; dense=true)`.
The default `dense=false` preserves the usual matrix storage.

The default `scalar=1` preserves frequency units. The device-aware
`sesolve(model, gate, initial, times)` handles this conversion, pulse stops, and
automatic energy centering. Returned states have their original lab-frame phase.

For delayed controls, `tstops=pulse_tstops(gate)` forces the adaptive solver to
resolve pulse boundaries and flattop transitions. Saving states on a dense time
grid alone does not force integration steps there. `get_unitary(model, gate)`
supplies these stops automatically. If supplying additional solver stop times,
combine them with `pulse_tstops(gate)`.

## Plot state trajectories

`state_amplitudes` computes complex overlaps with a dictionary of fixed reference
states. It accepts either a history of kets (or numerical vectors) or a
QuantumToolbox time-evolution solution. It preserves normalization and phase;
convert the amplitudes explicitly to populations with `abs2`.

```@example gates
using CairoMakie
references = Dict("Ground" => basis(2, 0), "Excited" => basis(2, 1))
amplitudes = state_amplitudes(references, result)
populations = Dict(key => abs2.(values) for (key, values) in amplitudes)
styles = Dict(
    "Ground" => (; color=:gray, linestyle=:dash),
    "Excited" => (; color=:purple, linewidth=2.5),
)
fig, ax = plot_trajectories(populations, result.times_states, styles;
    figure=(; size=(650, 350)),
    axis=(; xlabel="Time (ns)", ylabel="Population"),
    legend_options=(; position=:rt))
ylims!(ax, 0, 1.05)
fig
```

Use `result.times_states` for amplitudes: `result.times` can instead be the grid
used to sample expectation values. The value and style dictionaries must have
exactly matching keys, and every series must have the same length as the x vector.
Style entries accept ordinary Makie line options, as NamedTuples or Symbol-keyed
dictionaries; labels default to the dictionary keys converted to strings.

For a multi-panel figure, construct your axes and add curves to each one:

```@example gates
panels = Figure(size=(650, 600))
for (row, key) in enumerate(keys(populations))
    panel = Axis(panels[row, 1]; xlabel="Time (ns)", ylabel="Population", title=key)
    plot_trajectories!(panel, Dict(key => populations[key]), result.times_states,
        Dict(key => styles[key]); legend=false)
end
panels
```

The helpers also plot other real-valued trajectories, such as energies or pulse
coefficients. Precompute `real.(values)` or `imag.(values)` to plot amplitude
components. Reference states need not form a complete orthonormal basis, so their
squared overlaps need not sum to one.

## Overrides and endpoint checks

Gate parameters override matching model parameters. A parameter that changes an
idle coefficient must equal the idle value **exactly** at both endpoints; this is
an equality check, not a tolerance-based continuity test. An additional drive
parameter absent from the idle model has no idle-matching restriction.

Explicit gate duration must be finite, nonnegative, and contain every pulse;
`get_unitary` requires positive duration. Inferred timing is recomputed after
updates to pulse delay or duration. Scalar-only gates require explicit duration.
A zero additional Hamiltonian represents a gate that changes idle coefficients
only. `numerical` returns a static quantum operator when all coefficients are
static, and a `QobjEvo` when time dependence is present.

Named gates are stored in `model.gates[:name]`. A gate object is immutable;
use [parameter paths](parameters.md) to construct updated versions.

## Evolving several input states

Pass an ordered vector of kets to evolve their columns together. Automatic
energy centering subtracts the mean initial-state energy during integration,
then restores global phase at every saved time. Set `energy_shift=0` to disable
it. Native
QuantumToolbox calls retain angular-frequency units:

```julia
states = [model.states[(0, 0)], model.states[(1, 0)]]
solution = sesolve(numerical(model, gate; scalar=2pi), states, [0.0, gate.duration];
    progress_bar=false, tstops=pulse_tstops(gate))
M = get_gate_matrix(solution, states)
```

The evolved object is rectangular, with one column per input state. Projection
returns a small matrix in the chosen output basis. The three-argument `numerical` form supplies
`2π` and pulse event times internally:

```julia
M = numerical(model, gate, states; abstol=1e-8, reltol=1e-8)
fidelity = unitary_fidelity(sigmax(), M) # Includes relative phases by default.
probability_score = unitary_fidelity(sigmax(), M; include_phases=false)
leakage = 1 - sum(abs2, M.data) / length(states)
```

The two-argument `numerical(model, gate)` returns the Hamiltonian in frequency
units. Adding `states` returns a dimensionless `QuantumObject` operator after
full-model evolution and projection, retaining leakage. Its dimensions match the
selected basis; use `.data` to access the matrix entries.

Use `output_states` to choose another ordered orthonormal basis. Select
`frame=:interaction` explicitly to remove idle evolution before comparison.
`state_amplitudes` remains the helper for individual-ket trajectories.

Centering offsets use frequency units for device-aware evolution, gate matrices,
and calibration; the native-Hamiltonian batched `sesolve(H, states, times)`
overload uses angular-frequency units. `sesolveProblem` remains a native,
uncentered problem constructor. User integrator callbacks observe the centered
state; disable centering if a callback needs the original global phase.
