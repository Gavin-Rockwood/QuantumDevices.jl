# Model tomography

```@docs
Measurement
SequenceData
ModelTomographyProblem
predict_sequence
predict_sequences
reconstructed_model
tomography_residuals
```

## Data and losses

Scalar fitting uses `ModelTomographyProblem(model, paths; predict, observed,
uncertainty=1)`. `observed` must be a nonempty NamedTuple. `predict(model)` is
called once per evaluation and must return a NamedTuple with exactly the same
names; ordering may differ. Each value must be a finite real scalar. Residuals
follow the observation order.

`uncertainty` is a common positive scalar, or a NamedTuple with one positive
uncertainty for every observation. For example, use
`uncertainty=(ge=0.001, ef=0.002)` with `observed=(ge=4.8, ef=4.6)`.
The default of 1 gives unweighted least squares. Predictions, observations, and
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
norm equals the objective. Scalar fitting and sequence fitting use separate constructors.

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
