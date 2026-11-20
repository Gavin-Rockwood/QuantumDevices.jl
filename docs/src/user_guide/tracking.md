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
tracked = track_states(states; other=(energies=energies,))
@assert tracked.total_steps == length(values)
fig = Figure(size=(700, 320))
ax = Axis(fig[1, 1], xlabel="Detuning", ylabel="Tracked energy")
for label in 1:2
    lines!(ax, values, tracked.other.energies[label, :]; label="state $label")
end
axislegend(ax; position=:lt)
fig
```

Explicit `keys` can replace the default occupation tuples `(0,), (1,), ...`.
`save_steps=:end` returns only the final step, while still tracking through all
intermediate samples. Explicit step indices similarly select output, not execution.

## Dressed states

Recipes use `(;)` to form product states and turn on each node's
interaction along `trajectory`, default `x^3`, with 20 samples. Product keys such
as `(0, 1)` describe local energy levels; `model.states(0, 1)` gives a final
dressed ket and `model.energies(0, 1)` its energy. See [State labeling](state_labeling.md)
for ramp options, saved histories, transports, and diagnostics.

The default `track_states(...; method=:quick)` visits references in tracked order
and chooses the best still-available candidate, breaking ties by candidate order.
`method=:global` instead maximizes total squared overlap across the one-to-one
assignment. Use `(;method=:global)` to select it for a recipe ramp.

`accuracy.cumulative_log_overlaps` stores the logarithm of the product of matched
squared overlaps along the sampled path.
It is step-size dependent and is not a posterior probability or a fidelity to a
fixed reference state. Increase sampling density and check label stability, especially
near small gaps. Degenerate subspaces have no unique eigenvector labeling; individual
vectors can rotate even when the subspace itself varies smoothly.

Independent ramp diagonalizations can be batched per build with
`build(device, recipe; dressed_threads=4, parameters=Dict())`. The same keyword
is available when building a `GateRecipe` or `QuantumGate`, and when evolving a
`QuantumGate` with a device. Launch Julia with `--threads=4` (or more); workers
are capped by available default-pool threads and sample count.
The default `dressed_threads=nothing` honors each node's recipe
`dressing_options.thread_count` (default `1`). An explicit positive `Int`
overrides that count throughout the tree, without enabling disabled dressing or
changing saved recipes, drafts, captured definitions, or physical provenance.
Tracking inherits the effective dressed count at each node by default
(`tracking_threads=missing`). For example, `build(device, recipe; dressed_threads=4)`
uses up to four workers for both stages; `build(device, recipe; dressed_threads=4,
tracking_threads=2)` independently limits overlap work to two. Generic
`track_states(history; tracking_threads=2)` enables the same execution option.
Explicit tracking counts must be positive `Int` values.

Overlap matrices for adjacent column pairs are computed concurrently, followed
by the matrices bridging those pairs. Assignment, phase alignment, and diagnostics
are consumed in chronological order, preserving quick matching’s reference order
and global matching’s ordered behavior. Bounded windows retain at most twice the
worker count in reusable overlap matrices, plus the current selected-row matrix.
Each active worker owns reusable packing matrices for the two candidate collections
(and transported references when needed). Matrix multiplication writes overlaps
in place; buffers are owned by task slots rather than Julia thread IDs, so task
migration cannot share scratch storage. They do not retain dense overlaps for the
whole history. The serial quick path continues
to evaluate only available candidates. Parallel tracking evaluates dense pairwise
overlaps, which can cost more time and memory for small histories or many unused
candidates. Results are
numerically equivalent within tolerance; individual vectors near degeneracy
need not be identical. BLAS threading can compete with these tasks, so speedup
depends on the workload and BLAS settings. This setting changes neither Julia
nor BLAS thread settings.
