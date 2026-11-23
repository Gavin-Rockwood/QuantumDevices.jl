---
name: quantumdevices-tracking
description: Configure QuantumDevices recipe dressing and diagnose state labels across model parameter paths, basis changes, and sampled ramps.
---

Consult [state labeling](../../docs/src/user_guide/state_labeling.md) and [tracking](../../docs/src/user_guide/tracking.md). Inspect recipe component order, retained zero-based occupation tuples, dressing settings, and the physical ramp before changing labels. Distinguish occupation keys from array indices and energy ordering.

When recipe dressing is enabled, `dressing_options=(;)` selects the default product-state references and interaction-ramp settings at each node. Preserve existing dressing enablement and ramp settings; refine sampling and compare label stability near small gaps. Degenerate subspaces may be stable while individual eigenvectors rotate. Inspect model label diagnostics rather than treating cumulative overlap as fidelity.

```julia
# Assume device and recipe already define the intended model and dressing ramp.
model = build(device, recipe; dressed_threads=4, tracking_threads=2)
```

`dressed_threads` bounds independent ramp diagonalizations; `tracking_threads` bounds overlap work, each capped by available default-pool threads. Default `tracking_threads=missing` inherits the effective dressed count. These options do not change Julia or BLAS thread settings. Do not build merely to inspect the definition.

For existing chronological ket histories, use the focused [track_states skill](../quantumdevices-track-states/SKILL.md), which owns seed, subset, transport, auxiliary-data, result-shape, and assignment guidance. If models use changing construction coordinates, establish the correct basis transports before matching. For periodic-drive histories, use [Floquet basis and sweep](../quantumdevices-floquet-basis/SKILL.md).
