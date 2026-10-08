# Pulses, gates, and calibration

## Timed controls

```@docs
AbstractPulse
Pulse
pulse_tstops
DeviceGate
```

## Envelopes and ramps

```@docs
AbstractEnvelope
Constant
Gaussian
GaussianZero
GaussianSquare
Sech
Cosine
Blackman
Bump
ErfSquare
Slepian
DRAG
SineSquared
RampedFlattop
Envelope
envelope_value
validate_envelope
envelope_tstops
```

## Carriers

```@docs
AbstractCarrier
SineCarrier
IQCarrier
Carrier
carrier_value
```

## Evolution and metrics

```@docs
get_unitary
QuantumDevices.sesolve(::DeviceModel, ::DeviceGate, ::Union{QuantumDevices.QuantumObject,AbstractVector{<:QuantumDevices.QuantumObject}}, ::AbstractVector)
get_gate_matrix
state_amplitudes
plot_trajectories
plot_trajectories!
gate_unitary
unitary_fidelity
unitary_infidelity
gate_infidelity
```

## Calibration

```@docs
parameters
calibration_values
AbstractCalibrationSetup
CalibrationProblem
SciMLCalibrationSetup
calibration_problem
calibrated_gate
calibrate
```

`solve(problem, algorithm; kwargs...)` forwards directly to SciML's solver.
