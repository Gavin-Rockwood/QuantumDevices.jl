---
name: quantumdevices-spectral-tools
description: Choose a QuantumDevices spectral-analysis workflow from existing spectra and labels, sampled gaps, or a periodic-drive resonance or Stark-shift goal.
---

Start from the user's physical question and available data. Inspect the [models guide](../../docs/src/user_guide/models.md), [state-labeling guide](../../docs/src/user_guide/state_labeling.md), or [spectral reference](../../docs/src/reference/spectral_tools.md) only as needed.

- **Existing spectrum or state labels:** inspect the built model's `spectrum`, `energies`, and `states`, preserving component order and occupation labels. For labels across a parameter path, use the [tracking workflow](../quantumdevices-tracking/SKILL.md); for already collected ket histories, use [track_states](../quantumdevices-track-states/SKILL.md). Do not infer branch identity from energy order alone.
- **Already sampled crossing gaps:** use [fit_avoided_crossing](../quantumdevices-fit-avoided-crossing/SKILL.md). It requires a bracketed isolated minimum and at least five samples, without Floquet propagation.
- **Periodic-drive resonance:** use [find_resonance](../quantumdevices-find-resonance/SKILL.md) to track and fit two Floquet branches. Establish a genuinely periodic drive, reference states, frequency bracket, and units first. For modes or general parameter histories rather than a resonance fit, use [Floquet basis and sweep](../quantumdevices-floquet-basis/SKILL.md).
- **Drive-induced frequency shift:** apply `stark_shift` to a resonance result against its captured or explicit bare reference. It subtracts fitted drive frequencies; raw folded quasienergies are not shifts.

```julia
# Assume parameters and gaps are existing samples of an isolated crossing.
fit = fit_avoided_crossing(parameters, gaps)
fit.center, fit.minimum_gap, fit.rms_residual
```

Keep sample generation separate from analysis when samples already exist. Inspect continuity, competing modes, fit residuals, and interval dependence as appropriate; these do not certify gate fidelity or two-state isolation. Use [calibration](../quantumdevices-calibration/SKILL.md) when the actual goal is pulse optimization. Do not run model builds, scans, optimization, or evolution merely to inspect a setup.
