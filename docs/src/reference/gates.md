# Pulses, gates, and calibration

## Pulse interfaces

```@docs
AbstractPulse
AbstractParameterizedPulse
pulse_function
pulse_value
InternalPulseFunction
GenericPulseFunction
```

## Built-in pulses

```@docs
constant_pulse
gaussian_pulse
sine_pulse
sine_squared_pulse
ramped_flattop_pulse
available_ramps
```

## Gates and calibration

```@docs
DeviceGate
parameters
calibration_values
AbstractCalibrationSetup
SciMLCalibrationSetup
calibration_problem
calibrate
calibrated_gate
gate_unitary
unitary_infidelity
gate_infidelity
```
