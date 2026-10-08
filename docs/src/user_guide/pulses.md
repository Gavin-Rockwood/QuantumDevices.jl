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

## Quantum-control envelope library

All widths and ramp times below use physical time units. Amplitude and timing
stay on `Pulse`, so the same shapes work for drives and flux controls.

```@example pulses
shapes = (; gaussian=GaussianZero(0.15), gaussian_square=GaussianSquare(0.08),
    sech=Sech(0.15), cosine=Cosine(), blackman=Blackman(),
    erf_square=ErfSquare(0.08), slepian=Slepian(2.5; samples=129))
fig = Figure(size=(760, 380))
ax = Axis(fig[1, 1]; xlabel="Pulse-local time", ylabel="Dimensionless envelope")
for (name, shape) in pairs(shapes)
    lines!(ax, ts, [Pulse(shape; duration=1.0)(t) for t in ts]; label=String(name))
end
axislegend(ax; position=:rt)
fig
```

`GaussianZero` subtracts the centered Gaussian's endpoint baseline and normalizes
its peak to one. `GaussianSquare` uses these Gaussian edges around a flat plateau;
`ErfSquare` uses normalized error-function ramps. Both default to
`ramp_time=2sigma`, require both ramps to fit, and allow `ramp_time=0` for abrupt
edges. `Cosine` is a single cosine lobe, distinct from `SineSquared`.
`Blackman` is the standard symmetric Blackman window. These shapes have exact
zero endpoints when ramp times are positive. `Sech` and Slepian retain their
truncated endpoint values.

### Bump envelope

`Bump(k=2; center=nothing)` is a smooth compact-support envelope. With
`x = (t-center)/(duration/2)`, it equals `exp(k*x^2/(x^2-1))` inside
`abs(x) < 1` and zero outside. It peaks at one, and all derivatives vanish at
its support boundaries. Positive dimensionless `k` controls sharpness; larger
values narrow the peak.

```@example pulses
bump = Pulse(Bump(2); amplitude=0.2, duration=10.0)
bump(5.0) # unit envelope peak times amplitude
```

The default center follows half the pulse duration, so calibrating `duration`
automatically rescales the bump. Calibrate sharpness through `"drive/envelope/k"`.
An explicit `center` is in pulse-local time units; shifting it can truncate the
bump at the pulse window boundaries. `Bump()` also works as a reusable ramp in
`RampedFlattop(2.0; ramp=Bump())`.

### Slepian windows

`Slepian(time_bandwidth; samples=129)` uses the leading discrete prolate
spheroidal sequence (DPSS). The time-half-bandwidth product is dimensionless;
`samples` is an odd integer ≥3, with `0 < time_bandwidth < samples/2`.
The taper is peak-normalized and linearly interpolated over the pulse duration.
It is cached at construction; changing duration only rescales its clock.
Changing the time-bandwidth product rebuilds the cache.

```@example pulses
slepian_flux = Pulse(RampedFlattop(0.2; ramp=Slepian(2.5));
    amplitude=0.2, offset=0.05, duration=1.0)
@assert slepian_flux(0.0) == slepian_flux(1.0) == 0.05
adjusted = setpath(slepian_flux, "envelope/ramp_up/time_bandwidth", 3.0)
@assert getpath(adjusted, "envelope/ramp_up/time_bandwidth") == 3.0
```

Slepian endpoints are generally nonzero; `RampedFlattop` subtracts and normalizes
them when reusing the window as a ramp source. Only `time_bandwidth` is exposed
for calibration. Sample count is construction configuration; cached samples
are omitted from parameter discovery and save records.

### DRAG and I/Q modulation

`DRAG(sigma; beta=0)` produces the complex baseband envelope `G + im*beta*G′`,
where `G` is `GaussianZero` and `G′` is its analytic physical-time derivative.
`beta` has units of time and either sign is allowed. Normalization applies to
`G`, not the complex magnitude. The derivative quadrature can be nonzero at the
pulse boundaries even though `G` is zero there.

`IQCarrier(frequency; phase=0, reference=:pulse)` turns the complex envelope
`I + im*Q` into the real waveform `I*sin(θ) + Q*cos(θ)`. Its frequency uses cycles
per unit time, phase uses radians, and reference clock follows `SineCarrier`.
With zero Q it matches `SineCarrier`. IQ modulation requires real pulse amplitude
and offset; existing scalar-carrier and complex-baseband behavior is unchanged.

```@example pulses
drag_drive = Pulse(DRAG(0.15; beta=0.03); amplitude=0.02, duration=1.0,
    carrier=IQCarrier(5.0; phase=pi/2))
@assert drag_drive(0.3) isa Real
drag_gate = DeviceGate((; drive=drag_drive), param(:drive)*op(:q_charge))
changed_drag = setpath(drag_gate, "drive/envelope/beta", 0.04)
@assert getpath(changed_drag, "drive/envelope/beta") == 0.04
parameters(drag_gate)["drive/carrier/frequency"]
```

Use the real IQ-modulated coefficient with a Hermitian laboratory-frame drive
operator, such as charge. The shape alone is a complex baseband signal; it is
not automatically a Hermitian Hamiltonian. A DRAG envelope does not guarantee a
high-fidelity gate without calibrating the drive for the particular model.

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
adds these stops automatically, including flattop transitions. Interior joins and Slepian interpolation knots are included automatically.
Custom envelopes can implement `envelope_tstops(shape, duration)` to supply
pulse-local stop times. The device-aware `sesolve(model, gate, ...)` also adds
these stops automatically.

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

## Inspect controls in the REPL or a notebook

`display(pulse)` shows timing, amplitude, offset, envelope parameters, and carrier
information, followed by low-resolution UnicodePlots previews. The envelope plot
is dimensionless and excludes the carrier. Sine and IQ carriers have a separate constant
frequency plot in cycles per unit time; its phase is reported in radians. Custom
carriers show their parameters, since an instantaneous frequency cannot be inferred
from an arbitrary waveform.

```@example pulses
drive
```

`display(gate)` shows its Hamiltonian, scalar parameters, and each pulse on the
same gate-relative timeline, including delays and trailing idle time:

```@example pulses
gate
```

`print(pulse)` and `print(gate)` use one-line summaries. Compact contexts suppress
plots, and narrow displays retain the textual information. Complex baseband envelopes show I/Q labels with an IQ carrier.
These previews sample only the envelopes; they do not run a simulation.
