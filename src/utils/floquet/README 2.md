# Floquet utilities

Adapted from [SuperconductingCircuits.jl](https://github.com/Gavin-Rockwood/SuperconductingCircuits.jl/tree/main/src/Dynamics/Floquet).
All Hamiltonians and quasienergies use frequency units (cycles per unit time),
consistent with QuantumDevices. The evolution solver receives `2π * H`.

```julia
using QuantumDevices, QuantumToolbox

H = QobjEvo((0.17 * sigmaz(),
    (0.08 * sigmax(), (p, t) -> cos(2π * t))))
fb = get_floquet_basis(H, 1.0; abstol=1e-10, reltol=1e-10)
energies = fb.e_quasi
modes = fb.modes(0.25)
physical_state = exp(-2π * im * energies[1] * 0.25) * modes[1]

result = floquet_sweep(a -> QobjEvo((0.17 * sigmaz(),
        (a * sigmax(), (p, t) -> cos(2π * t)))),
    [0.0, 0.04, 0.08], 1.0;
    states_to_track=Dict(:zero => basis(2, 0)), use_logging=false)
tracked_energies = result["Tracking"].others["Quasienergies"]
```

`t0` chooses the absolute reference time for the basis. Mode sampling uses
`mod(t - t0, T)` and includes the quasienergy phase correction, making modes
periodic even for negative times. `H` must be periodic with period `T`.
Quasienergies are sorted within the principal zone; tracking reorders them by
state continuity but does not unwrap zone crossings.

Sweeps accept either one period or one period per point and reuse bases for
repeated parameter/period pairs. `sampling_times` defaults to `t0` and otherwise
must contain one time per point. `propagator_kwargs` forwards solver options;
direct solver keywords override those options. `tracking_kwargs` configures
the library's `track_states` function.

The qualified name `QuantumDevices.floquet_basis` remains available as an alias
of `FloquetBasis`. The low-level `propagate_floquet_modes(modes, H, t, T)` retains
the upstream within-period evolution behavior; pass `e_quasi` or use
`propagate_floquet_modes(fb, t)` for periodic Floquet modes.
