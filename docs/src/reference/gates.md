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
SineSquared
RampedFlattop
Envelope
envelope_value
validate_envelope
```

## Carriers

```@docs
AbstractCarrier
SineCarrier
Carrier
carrier_value
```

## Evolution and metrics

```@docs
get_unitary
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
