# Pulses, envelopes, and carriers

A `Pulse` owns its amplitude, baseline, duration, and delay. Evaluate it with
`pulse(t)`, where `t` is time since gate start. An envelope describes its shape;
a carrier supplies an oscillation when needed.

```@example pulses
using QuantumDevices, CairoMakie
flux = Pulse(RampedFlattop(0.2); amplitude=0.2, offset=0.05,
    duration=1.0, delay=0.1)
drive = Pulse(Gaussian(0.15); amplitude=0.02, duration=1.0,
    carrier=SineCarrier(5.0; phase=pi/2))
@assert flux(0.0) == flux(1.2) == 0.05
@assert flux(0.1) == flux(1.1) == 0.05
@assert flux(0.6) == 0.25
fig = Figure(size=(700, 340))
ax = Axis(fig[1, 1]; xlabel="Time", ylabel="Control value")
ts = range(0, 1.2; length=601)
lines!(ax, ts, flux.(ts); label="Flux")
lines!(ax, ts, drive.(ts); label="Drive")
axislegend(ax)
fig
```

During the inclusive active window `[delay, delay + duration]`, the value is
`offset + amplitude * envelope * carrier`. Without a carrier its multiplier is
one. Outside the window the value is `offset`. Duration must be positive and
delay nonnegative. Timing and Gaussian widths use physical time units.

`SineCarrier(frequency; phase=0)` uses cycles per unit time and radians. Its clock
starts at pulse onset. Choose `reference=:gate` to preserve a carrier's phase
reference to gate start when changing delay. Pulse offsets are not modulated.

## Reusable shapes

```@example pulses
shapes = (Constant(), Gaussian(0.15), SineSquared(), RampedFlattop(0.2))
[pulse(0.5) for pulse in (Pulse(shape; duration=1.0) for shape in shapes)]
```

`Constant()` is one throughout the window. `Gaussian(sigma; center=nothing)`
peaks at one; default center is half the duration. Its truncated endpoints are
not shifted to zero. `SineSquared()` gives one unit-height lobe with exactly zero
endpoints. `RampedFlattop(ramp_time)` has sin² rise and fall ramps by default.
Envelopes can also be inspected directly as `shape(local_time, duration)`.

## Independent ramps

Construct sources directly instead of passing constructors and keyword bundles:

```@example pulses
asymmetric = Pulse(RampedFlattop(0.2;
    rise_time=0.1, fall_time=0.3,
    ramp_up=Gaussian(0.2; center=0.4), split_up=0.4,
    ramp_down=SineSquared()); amplitude=0.2, duration=1.0)
@assert asymmetric(0.0) == asymmetric(1.0) == 0
@assert asymmetric(0.1) == asymmetric(0.7) == 0.2
parameters(asymmetric)["envelope/ramp_up/sigma"]
```

Named ramps are `:sine_squared`, `:linear`, and `:smoothstep`. The shared `ramp`
and `split` keywords provide defaults for `ramp_up`, `ramp_down`, `split_up`, and
`split_down`. Rise and fall times must fit within the pulse duration; zero time
means an abrupt edge on that side.

Envelope ramp sources use normalized source time `0:1`. Each side subtracts
its source endpoint and divides by the endpoint-to-split contrast. Source width
and center therefore use normalized units. The selected split must differ from
the relevant source endpoint. Normalization does not guarantee monotonicity or prevent
overshoot. Source parameters are stored explicitly and remain calibratable.

## Integration boundaries

Use `tstops=pulse_tstops(gate)` in direct `sesolve` calls. This prevents an adaptive
solver from stepping over a short delayed control. `get_unitary(model, gate)`
adds these stops automatically, including flattop transitions. Custom envelopes
with additional discontinuities need corresponding solver stop times.

## Gates and calibration

```@example pulses
gate = DeviceGate((; flux, drive), param(:drive) * op(:q_charge))
@assert gate.duration == 1.1
changed = setpath(gate, "drive/delay", 0.5)
@assert changed.duration == 1.5
@assert getpath(gate, "drive/carrier/frequency") == 5.0
parameters(gate)["drive/envelope/sigma"]
```

The gate infers its duration from the latest pulse end, including delay. This is
recomputed after pulse updates. An explicit gate duration adds trailing idle time
and constrains every pulse to fit. Flux controls use existing model parameter
names and must return to the model's idle value at gate endpoints. See
[gates](gates.md), [calibration](calibration.md), and [custom shapes](../development/extensions.md).
