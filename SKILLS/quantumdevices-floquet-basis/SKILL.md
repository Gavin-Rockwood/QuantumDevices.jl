---
name: quantumdevices-floquet-basis
description: Use get_floquet_basis, propagate_floquet_modes, propagate_floquet_states, and floquet_sweep for periodic-drive modes and histories.
---

Read the [Floquet guide](../../docs/src/user_guide/floquet.md), [basis implementation](../../src/utils/floquet/basis.jl), and [propagation](../../src/utils/floquet/propagation.jl).

`get_floquet_basis(H, period; t0=0, params=nothing, evolution_settings=(;), unitarity_atol=1e-7, mode_atol=1e-7, degeneracy_atol=1e-9, degeneracy_rtol=1e-10)` assumes a periodic drive in E/h units, a finite positive period, and absolute reference time `t0`. A finite pulse envelope is not automatically periodic. The solver supplies 2π once.

```julia
# Assume H is an existing periodic Hamiltonian in GHz and period is in ns.
basis = get_floquet_basis(H, period)
modes = propagate_floquet_modes(basis, H, period / 2)
physical_states = propagate_floquet_states(basis, H, period / 2)
```

`FloquetBasis` contains `e_quasi` in `[-1/(2period), 1/(2period))`, reference `states`, `period`, `t0`, `diagnostics`, and `settings`. Check unitarity, mode residual/orthogonality, and circular degeneracy diagnostics. Unacceptable defects throw; the implementation does not repair nonunitary propagation.

`propagate_floquet_modes(basis, H, t; params=nothing, evolution_settings=(;))` returns periodic modes at absolute time t. `propagate_floquet_states(basis, H, t; kwargs...)` adds the physical dynamical phase `exp(-2πim*e_quasi*(t-t0))`. Both require the same drive and parameters: basis stores data, not the drive closure.

For ordered parameter samples, [floquet_sweep](../../src/utils/floquet/sweep.jl) accepts `floquet_sweep(H_func, points, periods; t0=0, sample_times=nothing, params=nothing, seed=nothing, keys=nothing, save_steps=:all, floquet_threads=1, tracking_kwargs=(;), evolution_settings=(;), floquet_kwargs=(;))`. Factory runs once per point; scalar periods broadcast, otherwise supply one per point. Default sample times equal t0. It preserves point order and returns `FloquetSweep` with `points`, `periods`, `sample_times`, `tracking`, `diagnostics`, `settings`. A smaller seed selects branches; keys only label them. Threaded factories must own independent mutable workspaces; candidate histories are transiently retained. Use [track_states](../quantumdevices-track-states/SKILL.md) for tracking contracts and [find_resonance](../quantumdevices-find-resonance/SKILL.md) for a fitted resonance rather than a general sweep. Do not run propagation or scans merely to inspect setup.
