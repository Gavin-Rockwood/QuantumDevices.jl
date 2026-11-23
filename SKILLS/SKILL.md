---
name: quantumdevices
description: Navigate project-local QuantumDevices.jl skills for device construction, simulation, calibration, spectral analysis, tracking, and persistence.
---

This is the repository-local catalog, not an automatic global skill installation. Read only the smallest applicable skill below, then only references needed for the task. Existing user values, units, notebook conventions, solver settings, outputs, and ongoing edits take precedence over illustrative examples. These instructions do not authorize simulation by themselves; execute requested work, and do not launch builds, scans, evolution, or optimization for read-only inspection.

| Available artifact or intended result | Read |
| --- | --- |
| New/changed device, components, couplings, standalone recipe, model | [Models](quantumdevices-models/SKILL.md) |
| Pulse/control, prepared gate, time evolution, populations or gate matrix | [Controls](quantumdevices-controls/SKILL.md) |
| Built model or prepared gate; only compile H/runtime parameters | [make_hamiltonian](quantumdevices-make-hamiltonian/SKILL.md) |
| Existing evolution result | [Controls](quantumdevices-controls/SKILL.md), retrieval guidance only; skip evolution |
| Existing spectra/labels or uncertain spectral question | [Spectral routing](quantumdevices-spectral-tools/SKILL.md) |
| Existing sampled gaps | [fit_avoided_crossing](quantumdevices-fit-avoided-crossing/SKILL.md); skip sampling |
| Existing chronological ket histories | [track_states](quantumdevices-track-states/SKILL.md) |
| Recipe dressing, model labels or changing bases | [Tracking](quantumdevices-tracking/SKILL.md) |
| Periodic-drive resonance or signed Stark shift | [find_resonance](quantumdevices-find-resonance/SKILL.md) |
| Floquet modes or general periodic-drive sweeps | [Floquet](quantumdevices-floquet-basis/SKILL.md) |
| Fit pulse parameters/adopt calibrated gate | [Calibration](quantumdevices-calibration/SKILL.md) |
| Save/load devices, gates, attached artifacts | [Device I/O](quantumdevices-device-io/SKILL.md) |
| Documentation or package checks | [Docs/tests](quantumdevices-docs-tests/SKILL.md) |

Start from existing artifacts: a built model need not be rebuilt to inspect its energies, states or H; recorded results need no new dynamics. For implementation details use linked source/docs, not every skill. Reusable plumbing belongs in small library helpers; keep notebook cells focused on the physical workflow.
