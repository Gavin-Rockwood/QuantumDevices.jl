```@meta
CurrentModule = QuantumDevices
DocTestSetup = quote
    using QuantumDevices
end
```

# Envelopes

Envelope functions define the slowly varying amplitude profile used by drive coefficients.
Every built-in envelope follows the same calling convention:

```julia
envelope(t, drive_time; kwargs...)
```

The exported envelope functions are:

| Name | Registry key | Behavior |
|------|--------------|----------|
| `envelope_square` | `"square"` | Constant amplitude. |
| `envelope_gaussian` | `"gaussian"` | Gaussian pulse centered by default at `drive_time/2`. |
| `envelope_gaussian_ramp` | `"gaussian_ramp"` | Gaussian rise and fall with a flat top. |
| `envelope_sine_squared` | `"sine_squared"` | Sine-squared ramp segment. |
| `envelope_sine_squared_ramp` | `"sine_squared_ramp"` | Sine-squared rise and fall with a flat top. |
| `envelope_bump` | `"bump"` | Smooth compact bump. |
| `envelope_bump_ramp` | `"bump_ramp"` | Bump-shaped rise and fall with a flat top. |

## Registry-Based Construction

`get_envelope(drive_time, envelope_name, envelope_kwargs; digitize=false, step_length=2.3)` returns a one-argument function of time.
The `envelope_name` must be a key in `envelope_dict`.

```julia
using QuantumDevices

env = get_envelope(40.0, "sine_squared_ramp", Dict(:ramp_time => 10.0))
env(5.0)
```

Set `digitize=true` to sample the envelope on fixed steps before evaluation. The default `step_length` is `2.3`.

## API

```@docs
envelope_square
envelope_gaussian
envelope_gaussian_ramp
envelope_sine_squared
envelope_sine_squared_ramp
envelope_bump
envelope_bump_ramp
get_envelope
```
