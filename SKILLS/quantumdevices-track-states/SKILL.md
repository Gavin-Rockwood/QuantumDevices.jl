---
name: quantumdevices-track-states
description: Use track_states to match existing ket histories, preserve labels, reorder auxiliary data, and inspect continuity diagnostics.
---

Read the [tracking guide](../../docs/src/user_guide/tracking.md) and [implementation](../../src/utils/tracking/track_states.jl).

`track_states(state_history; seed=nothing, seed_transport=nothing, keys=nothing, other=(;), save_steps=:all, transports=nothing, method=:quick, degeneracy_atol=1e-9, degeneracy_rtol=1e-10, overlap_threshold=0.5, tracking_threads=1)` consumes ordered samples of normalized candidate kets. Coordinates must agree unless explicit QuantumObject transports map consecutive bases. Auxiliary NamedTuple histories have one aligned array per sample; `other.energies` enables degeneracy diagnostics.

```julia
# Assume state_history contains ordered candidate-ket samples and
# energy_history contains the corresponding candidate energies at each sample.
tracked = track_states(state_history; other=(energies=energy_history,), save_steps=:end)
final_energies = tracked.other.energies
```

The `TrackingResult` exposes `states`, `other`, zero-based candidate `indices`, `keys`, `saved_steps`, `total_steps`, `accuracy`, and full-path `diagnostics`. Multiple saved samples have key rows and saved-sample columns; one saved sample gives vectors. Ordinary indexing uses Julia array indices; callable access such as `tracked.states((0,))` uses keys (default `(0,), (1,), ...`).

A smaller `seed` selects a subset of candidates; keys alone do not. Seed is step zero outside the saved history; `seed_transport` requires seed and maps it to sample one. `save_steps=:end` still processes every sample. `:quick` greedily assigns in reference order; `:global` maximizes total squared overlap. Both phase-align matched kets. Inspect minimum overlaps, low-overlap flags, and degeneracy, and refine sampling to check stability. Cumulative overlap is step-size dependent, not fidelity; degenerate eigenvectors need not have unique labels. `tracking_threads` limits overlap workers, not Julia/BLAS configuration. For recipe-generated labels and model ramps, use the [tracking workflow](../quantumdevices-tracking/SKILL.md).
