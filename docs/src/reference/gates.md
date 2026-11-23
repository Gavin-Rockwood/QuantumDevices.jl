# Pulses, gates, and calibration

## Timed controls and envelopes

```@docs
AbstractPulse
AbstractPulseRecipe
PulseRecipe
PulseSequenceRecipe
Pulse
PulseSequence
pulse_tstops
pulse_duration
AbstractEnvelope
Envelope
EnvelopeSequence
Segment
InterpolatedEnvelope
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
QuantumDevices.envelope_value
QuantumDevices.validate_envelope
envelope_tstops
```

## Gate preparation and evolution

```@docs
Control
ControlInput
GateRecipe
PreparedGate
check_endpoints
GateEvolution
evolve
unitary
idle_phases
gate_matrix
```

## Calibration

```@docs
calibration_parameters
CalibrationProblem
CalibrationResult
calibrated_recipe
validate_calibration
QuantumGate
gate_infidelity
solve
```

`solve(problem; kwargs...)` uses the default bound-aware NLopt Nelder–Mead backend.
`solve(problem, algorithm; kwargs...)` selects a native SciML optimizer and
returns CalibrationResult; solving `problem.problem` returns the native solution.

## Static voltage control maps

```@docs
AbstractControlMap
AffineControlMap
TabulatedControlMap
map_values
check_control_map
add_control_map!
replace_control_map!
flux_map_from_frequency
control_map_diagnostics
```
