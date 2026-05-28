```@meta
CurrentModule = QuantumDevices
DocTestSetup = quote
    using QuantumDevices
    using QuantumToolbox
end
```

# Resonance and Calibration

QuantumDevices includes helper routines for finding driven resonances and calibrating pulse durations.
These tools are intended for small sweeps and calibration loops where the caller supplies the Hamiltonian, drive operator, reference states, and objective.

## Finding a Resonance

`find_resonance` has two entry points. The lower-level form accepts a function that maps each sampled parameter to a Hamiltonian:

```julia
result = find_resonance(H_func, freqs, reference_states)
```

The convenience form builds a sinusoidal drive around a static Hamiltonian and drive operator:

```julia
result = find_resonance(H0, drive_op, amplitude, freqs, reference_states)
```

`reference_states` is an array of states to track through the Floquet sweep. The return value is `[resonance_frequency, approximate_drive_time]`.

## Calibrating a Drive Duration

`calibrate_drive(drive_op, t_range0, psi0, to_min; kwargs...)` repeatedly samples a time interval and shrinks it around the best objective value.
The objective receives the final evolved state and returns a scalar to minimize.

```julia
to_min(psi_final) = 1 - abs(target' * psi_final)^2
best_time, best_value = calibrate_drive(drive_op, (0.0, 50.0), psi0, to_min)
```

The return value is `[best_time, best_value]`.

## FLZ Flat-Top Calibration

`get_FLZ_flattop` estimates a flat-top duration for a Floquet-Landau-Zener protocol from a static Hamiltonian, drive operator, frequency, amplitude, envelope function, ramp time, and two states.

This routine performs a Floquet sweep through the ramp and then searches phase samples to align the final state with the target Floquet superposition.

## API

```@docs
find_resonance
calibrate_drive
get_FLZ_flattop
```
