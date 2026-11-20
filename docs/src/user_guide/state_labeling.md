# Tracked data and dressed states

`model.eigensystem` is a `DressedStates` organization of the final eigenstates
and energies. Its immutable `.components` tuple gives the physical component
order used by occupation keys. For a target with order `(:q1, :coupler, :q2)`:

```julia
model.energies(1, 0, 0)       # occupation lookup, equivalent to ((1,0,0))
model.states(1, 0, 0)         # eigenvector in model construction coordinates
model.energies[1:5]           # ordinary one-based energy-ordered slice
model.eigensystem.components
```

Keys are plain zero-based integer tuples, not named labels. Component order
follows actual tensor inputs, recursively flattened: explicit grouped inputs
come first, then automatic inputs. It need not equal the original `components`
argument order. Complete keys are required. Selected energy levels retain their
original occupation keys rather than being renumbered.

`model.states`, `state(model, 1, 0, 0)`, and `states(model)` all return physical
eigenvectors in the model's construction basis, the same coordinates used by
stored `model.H`, `make_hamiltonian(model)`, and cached operators. The default static
`make_hamiltonian(model)` copies stored H; neither it nor the state API changes basis. At a leaf this is the operator
generator basis (the charge basis for built-in transmons). At a composite it is
the tensor product of explicitly selected upstream eigenvectors. Eigensolving
and labeling do not rotate the model's Hamiltonian or operators into a new basis.

`stateindex(model, 1, 0, 0)` returns a zero-based **energy-order index**, not a
coordinate index; `stateenergy(model, 1, 0, 0)` is the corresponding energy lookup.
Thus `model.H * state(model, key)` equals `stateenergy(model, key) * state(model, key)`
within numerical tolerance. Do not replace these eigenvectors with coordinate
`basis` kets unless the model construction Hamiltonian is diagonal.

An explicit pure reduction selects upstream dressed eigenvectors as its new
basis. Its Hamiltonian is naturally diagonal and its exposed eigenstates are
coordinate basis kets in ascending energy order, with original labels and
inherited tracking diagnostics. No new tracking ramp or eigensolve occurs there.
Other builds still eagerly provide eigensystems to supply spectra, labels, and
subsequent state selections; there is no general lazy eigensystem interface.
Known diagonal product stages reuse their energies and basis kets without another
eigensolve. Basis projectors are transient; no duplicate dense basis map is retained.

## Array containers with shared tracking

States, energies, and tracking-quality arrays are immutable `TrackedData`
containers. Their array and QuantumObject contents have ordinary mutable behavior;
there is no deep immutability promise or custom mutation policing. Normal Julia
indexing, including linear matrix indexing, is one-based. Ordinary slices, `copy`,
`collect`, broadcasts, and arithmetic yield editable arrays. The wrapper provides
no assignment forwarding or mutation guarantees; views use normal AbstractArray
machinery rather than a special copying rule.

Every wrapper points to the same owner through `.tracking`:

```julia
model.states.tracking === model.energies.tracking
model.eigensystem.tracking.keys
model.eigensystem.tracking.accuracy.minimum_overlaps(1, 0, 0)
model.eigensystem.tracking.diagnostics(1, 0, 0)
```

`model.spectrum` remains an energy-array compatibility alias. `model.labeling`
is a migration alias for the tracking owner when available; new code should use
`model.eigensystem`. For `dressing=false`, states and energies still support
normal indexing, but occupation lookup raises `ArgumentError`; keys are empty,
and accuracy/final diagnostics are `nothing`.

## Generic tracking and saved samples

```julia
result = track_states(state_history;
    keys=(:left, :right),
    other=(energies=energy_history, tags=tag_history),
    save_steps=:all)
result.states(:left)
result.other.energies(:left)
result.other.tags(:left)
```

Each input sample contains normalized QuantumToolbox kets. Subsequent samples
may have extra candidates. `other` is a NamedTuple of candidate-indexed histories;
every field follows exactly the state assignment and preserves its element type.
Only the optional `energies` field has special validation: finite real values
supply near-degeneracy diagnostics. Without energies, only overlap diagnostics
are available. Default keys are `(0,), (1,), ...`, so `result.states(0)` means
`result.states((0,))`. Explicit hashable keys such as Symbols also work. Keys must be unique, and callers
are responsible for keeping their hash/equality behavior stable while in use.

`save_steps` accepts `:all`, `:end`, or nonempty strictly increasing one-based
sample indices. Every input step is processed, including unsaved intermediates.
With multiple saved samples, arrays have key rows and saved-sample columns. With
exactly one saved sample, all saved fields are vectors: callable lookup returns
a scalar element instead of a row history. The owner stores immutable `.keys`
and `.saved_steps`, together with `.total_steps`.

```julia
final = track_states(state_history; other=(energies=energy_history,), save_steps=:end)
final.states(0)                    # one QuantumObject
final.other.energies(0)            # one energy
final.accuracy.minimum_overlaps(0)
final.saved_steps
```

Accuracy fields are `overlaps`, `cumulative_log_overlaps`, `minimum_overlaps`,
`low_overlap`, `degenerate`, and `ambiguous`. Cumulative fields include unsaved
intermediate steps. `.diagnostics` summarizes the entire processed path even
when the last saved sample is earlier than the final input sample. Overlaps are
squared magnitudes, and diagnostics are continuity measures, not calibrated
probabilities that a label is correct.

Histories share coordinates unless supplied with one QuantumObject linear
`transport` per transition through `transports=...`. Each map takes previous
coordinates into current coordinates. Projection norm loss is retained, not
normalized away. The default `method=:quick` visits references in tracked order and selects each
reference's best still-available candidate, with the first candidate winning ties.
Set `method=:global` to maximize total squared overlap under a one-to-one
assignment. Both methods phase-align selected states without mutating the input. Auxiliary
values are reordered, never phase-transformed. There is no automatic transport
between independently rebuilt models: each build establishes its own ramp keys.

### Following a small seed through larger candidate samples

Pass `seed` to supply a separate reference collection before the first history
sample. For example, follow four reference kets through samples containing 125
candidate eigenstates in the same coordinate basis:

```julia
using QuantumToolbox: basis
seed = [basis(125, i) for i in 0:3]
result = track_states(candidate_history;
    seed,
    keys=(:reference0, :reference1, :reference2, :reference3),
    other=(energies=energy_history,),
    save_steps=:all,
    tracking_threads=4)
result.states(:reference0)
result.indices(:reference0)  # zero-based candidate indices at saved samples
```

Each seed ket must be normalized. Keys must be unique and match the seed count;
omitting keys gives `(0,), (1,), ...` for the seed trajectories. Each candidate
sample must contain at least as many kets as the seed, but candidate counts may
vary between samples. Samples and candidate-indexed auxiliary histories are
unchanged: matching retains only the selected trajectories in the result.

The seed is a step-zero reference outside the supplied history. It is matched to
the first sample, and subsequent matches use the previously selected,
phase-aligned kets. `save_steps=[1]` saves the first supplied sample, and
`total_steps` counts supplied samples, excluding the seed. Auxiliary histories
need no seed entry. First-sample overlap accuracy and cumulative diagnostics
include the seed-to-first match, including any projection norm loss.

If seed coordinates differ from the first sample, supply a QuantumObject linear
map as `seed_transport=...`. Without it, their bases must agree. This map is
separate from `transports`, which still contains one map between each pair of
successive history samples. Supplying `seed_transport` without a seed is an error.
Neither transport renormalizes projected reference kets.

Seeded transitions run sequentially because their references depend on the
previous assignment. Within each transition, `tracking_threads` splits candidate
columns across default-pool workers with independent reusable scratch. Four
references against 125 candidates therefore compute a 4-by-125 overlap matrix,
not a 125-by-125 matrix. Both matching methods retain their existing deterministic
assignment order and phase alignment. With no seed, the existing adjacent-sample
parallel strategy and first-sample identity assignment are unchanged.

## Interaction-ramp meaning and ambiguity

Leaves assign singleton occupation tuples to their own energy levels. At each
combine, product keys follow only that node's new interaction along `H0 + g V`;
upstream interactions remain incorporated in child eigenstates. By default,
`dressing_options=(;)` uses `method=:quick`, 20 samples including endpoints, `g=x^3`, degeneracy
`tolerances atol=1e-9, rtol=1e-10`, and an overlap warning threshold of `0.5`.

```julia
options = (;method=:global, nsteps=60, trajectory=identity, overlap_threshold=0.7)
recipe = Recipe(device, components; dressing=true, dressing_options=options)
unlabeled = Recipe(device, components; dressing=false, dressing_options=options)
```

Trajectories require finite nondecreasing real samples in `[0,1]`, with exact
endpoints zero and one. Zero total new interaction uses only one final matching
step, using the known diagonal product energies and basis kets without another
eigensolve. A labeled parent requires labeled upstream composites. A key follows the
sampled ramp rather than the largest final bare character. Increase `nsteps` and
inspect stability where it matters.

Degeneracies satisfying `abs(Ei-Ej) <= atol + rtol*max(abs(Ei),abs(Ej))` are flagged.
Deterministic assignment and phase alignment do not make individual labels unique
inside a degenerate subspace. No subspace mixing or adaptive retries are applied.
Inherited minimum overlaps, cumulative logs, and ambiguity propagate through the
recipe tree; cumulative logs include that preparation history and depend on the
tree as well as the path. Final dressed rows, keys, auxiliary data, and diagnostics
are all organized in ascending energy order, without changing generic tracker
key order or recipe energy-index selections.

Independent ramp diagonalizations can be batched per build with
`build(device, recipe; dressed_threads=4, parameters=Dict())`. The same keyword
is available when building a `GateRecipe` or `QuantumGate`, and when evolving a
`QuantumGate` with a device. Launch Julia with `--threads=4` (or more); workers
are capped by available default-pool threads and sample count.
The default `dressed_threads=nothing` honors each node's recipe
`dressing_options.thread_count` (default `1`). An explicit positive `Int`
overrides that count throughout the tree, without enabling disabled dressing or
changing saved recipes, drafts, captured definitions, or physical provenance.
Samples retain ramp order and state tracking remains sequential. Results are
numerically equivalent within tolerance; individual vectors near degeneracy
need not be identical. BLAS threading can compete with these tasks, so speedup
depends on the workload and BLAS settings. This setting changes neither Julia
nor BLAS thread settings.

`dressing` must be a `Bool`. `dressing_options` must be a named tuple;
it is preserved in snapshots, drafts, and saved recipes even when `dressing=false`,
when it is unused during model construction. The recipe keywords replace the former
`labeling` keyword; `model.labeling` remains the tracking diagnostics property.

Unknown option names and invalid scalar values are rejected even when dressing is
disabled. Trajectory samples are generated and validated during an enabled build;
no derived sample array is stored on the recipe. Options merge over the defaults,
so omitted keys retain their default values. Callables retain their identity and
may capture mutable external state outside the snapshot boundary.
