---
name: quantumdevices-find-resonance
description: Use find_resonance for Floquet-backed periodic-drive resonance estimates and stark_shift for their signed frequency shifts.
---

Read the [Floquet guide](../../docs/src/user_guide/floquet.md), [resonance implementation](../../src/utils/spectral_tools/resonance.jl), and [model adapter](../../src/utils/spectral_tools/model_adapters.jl).

Choose the available input form:

- `find_resonance(H_func, frequencies, reference_states; keys=nothing, t0=0, params=nothing, photon_order=1, reference_frequency=nothing, floquet_threads=1, tracking_kwargs=(;), evolution_settings=(;), floquet_kwargs=(;), fit_kwargs=(;))`: factory `H_func(f)` returns a drive with fundamental period `1/f`; exactly two normalized orthogonal reference kets select branches.
- `find_resonance(H0::QuantumObject, drive::QuantumObject, amplitude, frequencies, reference_states; phase=0, kwargs...)`: constructs `H0 + amplitude*cos(2π*f*t+phase)*drive`; phase uses absolute time. Operators must be compatible and Hermitian.
- `find_resonance(model::NumericalModel, drive::Sym, amplitude, frequencies; states, phase=0, photon_order=1, kwargs...)`: two distinct complete occupation tuples supply model reference states and a bare reference frequency `abs(ΔE)/photon_order`. The symbolic drive must have cached complete operator products in the model's scope; it does not rebuild.

```julia
# Assume model is built, drive_operator is its compatible symbolic drive,
# and frequencies brackets an isolated crossing with at least five samples.
result = find_resonance(model, drive_operator, amplitude, frequencies;
    states=[(0,), (1,)])
shift = stark_shift(result)
```

Frequencies are positive finite distinct reals, sorted internally. H uses E/h units; propagation supplies 2π once. The method tracks two Floquet branches, fits their circular quasienergy gap, and requires a resolved bracketed crossing. `photon_order` is a positive integer describing the reference transition; the drive period remains `1/f`.

`ResonanceResult` owns `sweep`, `fit`, `reference_frequency`, `photon_order`, `metadata`, and `diagnostics`; convenience properties include `frequency`, `minimum_gap`, `drive_time=1/(2minimum_gap)`, `gaps`, `tracking`, `state_keys`, `quasienergies`, and `frequencies`. Inspect fit residuals, overlaps, degeneracy, and potential competing modes. The time is a two-state estimate, not a calibrated pulse duration or fidelity certificate.

[stark_shift](../../src/utils/spectral_tools/stark_shift.jl) returns `result.frequency-reference_frequency`, without solving again. Low-level results need `stark_shift(result; reference_frequency=bare_frequency)` or `stark_shift(result, model; states=result.state_keys, photon_order=result.photon_order)`. Folded quasienergy differences alone are not Stark shifts. If gaps are already sampled, use [fit_avoided_crossing](../quantumdevices-fit-avoided-crossing/SKILL.md) directly. Do not launch scans just to inspect definitions.
