# Model tomography

```@docs
Measurement
SequenceData
ModelTomographyProblem
predict_sequence
predict_sequences
reconstructed_model
tomography_residuals
physical_parameters
```

## Data and losses

Sequence prediction uses plain sequences and shared settings:

```julia
observed = predict_sequences(target, sequences; preparation, measurement)
problem = ModelTomographyProblem(initial, paths;
    sequences, preparation, measurement, observed, gates=initial.gates)
```

`observed` is one probability vector per sequence, ordered by measurement outcome.
For integer counts, use `counts=measured_counts` instead; do not provide both.
`shots` can be a shared positive integer or a vector matching the sequences.
Vector entries may be `nothing`; mixed shot qualification requires explicit
`loss=:least_squares`. Counts infer totals and check any supplied shots.
`predict_sequence(model, sequence; preparation, measurement)` handles one sequence.
Prediction requires no observed distributions. Retain `SequenceData` and the
three-argument problem constructor for per-experiment preparation or measurement.
Do not mix scalar `predict` with `sequences`.

Scalar fitting uses `ModelTomographyProblem(model, paths; predict, observed,
uncertainty=1)`. `observed` must be a nonempty NamedTuple. `predict(model)` is
called once per evaluation and must return a NamedTuple with exactly the same
names; ordering may differ. Each value must be a finite real scalar. Residuals
follow the observation order.

`uncertainty` is a common nonnegative scalar, or a NamedTuple with one nonnegative
uncertainty for every observation. For example, use
`uncertainty=(ge=0.001, ef=0.002)` with `observed=(ge=4.8, ef=4.6)`.
Zero uncertainty leaves that observation unweighted; it does not impose an exact
equality constraint. Thus `uncertainty=0.0` gives unweighted least squares, like
the default of 1. Predictions, observations, and
uncertainties must use the same units; errors are independent.
`loss=:auto` and `:least_squares` both minimize uncertainty-scaled residuals.

For `SequenceData`, counts infer their shot total. Probabilities may include
`shots` without rounding into reconstructed counts. Keep observations in the
measurement's outcome order; labels are not array indices.
`loss=:auto` uses multinomial negative log likelihood when every record carries
shots, and unweighted least squares when none does. Mixed shot-qualified records
require explicit `loss=:least_squares`; shot totals do not weight that loss.
Likelihood omits parameter-independent constants and uses `probability_floor`
only inside logarithms. Zero-count terms contribute zero.

`tomography_residuals(setup, values)` is available for least squares; its squared
norm equals the objective. Scalar and direct sequence fitting share the two-argument constructor.

## Measurements and prediction

Use `Measurement(effects; labels=..., acts_on=...)` for a complete POVM or
`Measurement(observable; acts_on=...)` for eigenspace projectors. Degenerate
observable eigenstates share an outcome. Effects must sum to identity; include
leakage as additional outcomes. Subsystem indices are one-based and ascending.
Preparations can be kets, density matrices, or functions of the candidate model.
Measurements can also be functions returning `Measurement`.

Fixed states and effects are interpreted in the candidate's retained component
basis. They are not automatically transformed from an external laboratory basis.
Use model-dependent functions when physical updates require regenerating them.

`predict_sequence` needs no observed data. `predict_sequences` ignores observed
distributions and shares each required gate propagator across all records in that
call. Gate clocks reset at each gate's local origin; empty sequences measure the
preparation directly. No propagator cache survives between candidates.

## Solver options

Bounds such as `lb` and `ub` go to `ModelTomographyProblem`; select an optimizer
that supports them. Optimizer settings and callbacks go to `solve`. The wrapper
returns the native SciML solution. Automatic differentiation through physical
reconstruction is not guaranteed; derivative-free optimization is a useful start.

`threaded=true` builds gate propagators concurrently. Shared controls and solver
callbacks must be safe for concurrent calls. `evolution_kwargs` forwards native
ODE solver settings. Both options affect sequence fits only.
Selected slash paths must refer to distinct writable scalar inputs; aliases and
overlapping paths are rejected. Inspect solver retcodes, held-out predictions,
other initial guesses, and dressed-state tracking confidence before accepting a
fit. Low loss does not establish unique parameter recovery.

## Projected sequence fitting

Supply `subspace` as an ordered orthonormal ket vector or a function returning
one for the candidate model, for example:

```julia
subspace = m -> [m.states[label] for label in ((0,0), (1,0), (2,0))]
setup = ModelTomographyProblem(model, data, paths;
    subspace, leakage_outcome=4, threaded=true)
```

Only the selected columns are integrated in the full Hamiltonian. The projected
gate matrices are cached once per candidate and multiplied to evaluate sequences.
The preparation must lie in the supplied basis. Projection after each gate loses
paths that leave the basis and later return; compare against full-space predictions
and enlarge the basis when needed. This option applies to sequence fits only.

`leakage_outcome` is a one-based index into each measurement's ordered effects.
That effect must be the complement of the subspace projector.
Discarded probability is assigned to this absorbing outcome without conditioning
on survival. Without it, appreciable discarded probability raises an error.
`predict_sequence` and `predict_sequences` accept the same keywords. The default
`subspace=nothing` retains full-space evolution.

## Parameter normalization

Fits default to dimensionless optimizer coordinates, with
`physical = initial + scale .* coordinates`. Initial coordinates are zero;
scales are absolute initial parameter values, or 1 for zero-valued parameters.
Use `parameter_scales=[...]` in path order to set characteristic changes, or
`scale_function=(path, initial) -> ...`. Scales must be finite and positive and
are fixed for the entire fit. Choose an explicit scale for small shifts starting
at zero. Parameter scaling does not change observation uncertainty weights.

`solution.u`, residual inputs, and reconstruction inputs are optimizer
coordinates. Use `physical_parameters(setup, solution.u)` to report physical
values. Bounds `lb` and `ub` are supplied in physical units and transformed for
the optimizer. Set `normalize=false` to use physical coordinates directly.
Optimizer simplex sizes and coordinate tolerances use the chosen coordinates;
objective tolerances such as `g_abstol` keep their existing meaning.

Component reconstruction uses `threaded_components=true` by default in both
tomography constructors, independently of `threaded` gate propagation. Changed
component constructors run concurrently when multiple components change and
Julia has multiple threads. The normal builder assembles operators and tracks
dressed states using the requested reference and intermediate steps. Use
`threaded_components=false` for serial reconstruction or custom constructors
that are not thread-safe. Start Julia with `--threads=auto` to enable workers.

With `normalize_probabilities=false`, full-space prediction corrects small state-normalization drift (at most `1e-6`)
after each gate, preventing accumulation when a numerical propagator is reused
many times. Larger per-gate errors raise an error asking for tighter evolution
tolerances. Projected-subspace leakage is never normalized away.

Use `normalize_probabilities=true` in prediction or sequence tomography to normalize
finite, nonnegative predicted outcome weights. Projected missing population is
assigned to the complementary leakage outcome first, preserving physical leakage.
Projected loss without a leakage outcome still raises an error, as do negative
weights beyond rounding tolerance, non-finite weights, or zero total weight.
The default `false` retains strict normalization checks. This differs from
`normalize`, which scales optimization parameters. Probability normalization
does not correct integration errors in relative phases or outcome proportions.

## Reconstruction and state references

Every physical candidate uses `make_model`, without an incremental reconstruction
cache. Tomography defaults to `reference_model=:initial, tracking_steps=0`: one
candidate diagonalization, matched to the fixed initial eigenbasis. Transmon
basis maps carry overlaps into the same physical charge basis after EC/EJ changes.
Reference states are not frozen and energies always come from the candidate.
`tracking_steps=N` adds N intermediate models along numeric parameter changes;
model structure must stay fixed. Confidence is accumulated squared physical
overlap along this path; it is not a correctness probability. Inspect it and
increase steps near avoided crossings. Candidate evaluation order has no effect.
Automatic initial references fall back to bare tracking when basis transport is
unavailable; incompatible explicit references raise an error.

Use `reference_model=nothing, tracking_steps=18` for the original bare tracking.
Normal `make_model` builds default to that bare trajectory, and accept
`reference_model=initial, tracking_steps=0` for a direct reference match.

Physical rebuilding does not inherit gates. Pass `gates=baseline.gates` explicitly
to sequence fitting and, when needed, prediction. Gate-only updates reuse idle
physics; selected gate-parameter paths are explicit updates and remain supported.
Use `Component(component; operators=...)` to extend operators while preserving
its derived basis metadata. Six-argument custom construction remains supported.

### Timing sequence fits

Pass `diagnostics=report -> ...` to `ModelTomographyProblem` or
`predict_sequences` to receive timings in seconds after each successful prediction:

```julia
reports = []
# Add this keyword to your sequence ModelTomographyProblem call:
# diagnostics = report -> push!(reports, report)
# After evaluating or solving:
# reports[end]
```

A sequence fit report contains `rebuild_seconds`, `gate_update_seconds`,
`setup_seconds` (sequence inputs and subspace), `propagation_seconds` (all gate
builds), `gate_seconds` (gate-name => seconds pairs), `sequence_seconds`
(sequence multiplication and measurements), `total_seconds`, and `parameters`
(in physical units). Direct `predict_sequences` reports omit rebuilding,
gate updates, and parameters. The total excludes the reporting callback and
loss reduction. Initial constructor validation also produces a report.
First-call timings include compilation; compare subsequent evaluations.
In threaded runs, individual gate times overlap: compare their wall-clock
`propagation_seconds` rather than summing `gate_seconds`. Reports arrive on the
calling task after gate workers finish. Diagnostics are disabled by default;
callback errors propagate to the caller.

### Evolution storage

Model/gate evolution uses dense Hamiltonian and initial-state matrices by default.
Use `evolution_kwargs=(; dense=false)` to preserve their original storage.
The same options apply to `numerical(model, gate)`, `get_unitary(model, gate)`,
`numerical(model, gate, states)`, and `sesolve(model, gate, initial, times)`.

For a supplied Hamiltonian, `get_unitary(H, duration; dense=true)` also converts
static matrices or the matrix terms of a `QobjEvo` to dense storage before solving.
Time-dependent coefficients are preserved, and the input Hamiltonian is not modified.
Use `dense=false` to preserve its original storage.
