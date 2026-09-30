# Tracking states and dressed labels

Energy ordering alone changes at crossings. `track_states` associates states
between successive parameter values by squared overlaps, retaining initial labels.
Auxiliary histories such as energies are reordered using the same matches.

```@example tracking
using QuantumDevices, QuantumToolbox, CairoMakie
values = range(-1, 1; length=41)
spectra = [eigenstates(0.5δ * sigmaz() + 0.1 * sigmax()) for δ in values]
energies = [s.values for s in spectra]
states = [begin
    λ, ψ = s
    ψ
end for s in spectra]
tracked = track_states(states; other_sorts=Dict("energy" => energies))
@assert length(tracked.confidence) == length(values)
fig = Figure(size=(700, 320))
ax = Axis(fig[1, 1], xlabel="Detuning", ylabel="Tracked energy")
for label in 1:2
    lines!(ax, values, [step[label] for step in tracked.others["energy"]]; label="state $label")
end
axislegend(ax)
fig
```

A labeled dictionary can replace the default first-step index labels.
`return_steps=:end` returns only the final step, while still tracking through all
intermediate samples. Explicit step indices similarly select output, not execution.

## Dressed states

`get_dressed_states(H0s, Hint)` starts with local energy eigenstates, forms their
product states, and turns on the interaction along `trajectory`, default `x^3`.
`make_model` performs this automatically. Product labels such as `(0, 1)` describe
local energy levels; `model.states[label]` gives a final dressed ket and
`model.others[label]` its energy.

The default `QuantumDevices.expansive_tracking` repeatedly chooses the largest
available overlap across unmatched pairs. `QuantumDevices.quick_tracking` assigns
in initial-state order. Both are greedy, not globally optimal assignment solvers.
Custom methods return `(indices, squared_overlaps)` in the input-state order.

Confidence is the product of matched squared overlaps along the sampled path.
It is step-size dependent and is not a posterior probability or a fidelity to a
fixed reference state. Increase sampling density and check label stability, especially
near small gaps. Degenerate subspaces have no unique eigenvector labeling; individual
vectors can rotate even when the subspace itself varies smoothly.
