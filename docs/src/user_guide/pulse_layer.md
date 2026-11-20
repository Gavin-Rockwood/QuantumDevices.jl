# Pulses and normalized envelopes

Standard shape constructors return a common `Envelope` with a physical reference
`duration`. Constructor widths, centers (where supported), interpolation nodes and
ramp times use that same physical time unit. The shape payload stores normalized
coordinates on `[0,1]`; amplitudes may be negative or complex.

```julia
using QuantumDevices
e = Gaussian(; duration=20., sigma=4.)
pulse = Pulse(e; amplitude=.2, offset=.05, delay=3.)
stretched = Pulse(e; duration=40.) # effective sigma is now 8
```

Gaussian is centered at half its window. There is no Gaussian `center` keyword:
`Pulse.delay` places the whole waveform. Sech and Bump retain physical `center`
keywords. `GaussianZero`, `GaussianSquare`, `ErfSquare`, `SineSquared`, `Cosine`,
`Blackman`, `Slepian`, `DRAG` and `RampedFlattop` also require a duration.
Bump takes physical support half-width and optional dimensionless `k`.

Duration-only edits stretch the entire normalized waveform, including widths,
interpolation nodes and DRAG beta. `Pulse` forwards its envelope's duration;
there is no second Pulse duration field. Explicit overrides create a new wrapper
sharing its shape payload and leave the original unchanged. To change duration
while holding physical sigma or beta fixed, supply both fields in one calibration
candidate: all final physical fields are normalized against the final duration.

## Regions and short compositions

`e.rise` / `e.fall` contain normalized interval metadata or `nothing`.
`e.rising` / `e.falling` select those actual regions, sharing the original payload.
The default view duration is the selected fraction of the reference duration.
Gaussian, SineSquared, Blackman and Slepian provide geometric half-window views;
these are not a guarantee of monotonicity. Missing region requests throw.
Public constructor intervals use physical times and are normalized once.

```julia
window = SineSquared(; duration=10.)
e = EnvelopeSequence(
    (shape=window.rising,),
    (shape=Constant(; duration=20.),),
    (shape=window.falling,),
    (shape=Constant(; duration=2.), scale=0),
)
pulse = Pulse(e; amplitude=.2)
```

Pieces accept `shape`, optional `duration`, and optional finite real or complex
`scale`. `Segment(shape; duration=nothing, scale=1)` is the equivalent record.
A different piece duration stretches an ordinary envelope or region.
EnvelopeSequence duration is the sum of its pieces: whole-sequence stretching is
rejected, including for a nested sequence piece. Internal joins belong to the
next piece; the final piece owns the final endpoint. Piece delays, carriers and
offsets are deliberately absent.

`RampedFlattop(5.; duration=30.)` evaluates rise, plateau and fall directly.
Built-in ramps are `:linear`, `:smoothstep`, and `:sine_squared`. Custom source
regions are sampled forward on their unit domain; use `.rising` and `.falling`
explicitly rather than relying on automatic reversal of a source shape.

## Independent pulses and clocks

`PulseSequence(pulses...; offset=0, delay=0)` combines independent full controls.
The first child delay is measured from sequence onset; subsequent delays are
gaps after the previous child. Active children supply their full offset; gaps
return the sequence offset. The enclosing sequence duration includes child gaps
and excludes its own delay. Nested clocks remain local to the immediate parent.

`Pulse(e, frequency; phase=0, reference=:pulse, ...)` mixes complex I+iQ as
`offset + amplitude*(I*cos(θ)-Q*sin(θ))`; amplitude and offset must be real.
Frequency is cycles per time unit and phase is radians. `reference=:sequence`
uses the immediate parent's clock, including the child's delay; standalone it
uses input time. Outside a Pulse's closed window, evaluation returns its offset.
`pulse_tstops` returns boundaries, joins and knots shifted by placement.

## Custom functions, parameters and saved pulses

```julia
custom = Envelope((p,u)->p.level*u^2, (;level=.4); duration=10., rise=(0.,10.))
changed = setpath(custom, "parameters/level", .8)
```

Custom callbacks receive `(parameters, normalized_time)`. Constructor endpoint
validation can call the callback; evaluation checks finite scalar output.
There is no callable-arity guessing. To migrate an old physical callback
`oldf(p,t,d)`, explicitly use `(p,u)->oldf(p,u*T,T)` with `duration=T`.
DRAG uses normalized derivative `g′(u)` with coefficient `beta/T`; stretching
preserves the complete I/Q field. Region selection samples that field and does
not insert an extra derivative factor.

`parameters` and `calibration_parameters` expose physical slash paths such as
`envelope/sigma`, `envelope/times/2`, and `envelope/segments/1/scale`; caches and
normalized shape fields are excluded. Pulse `duration` and `envelope/duration`
refer to one timing target; duplicate calibration writes through aliases fail.
Derived sequence totals cannot be selected independently. Current calibration
coordinates remain finite real scalars; complex scales can be edited directly.

`save_pulse` / `load_pulse` use versioned records. New wrapper, shape and region
IDs store normalized configuration and reconstruct caches. Old v1 physical-time
records load into private compatibility envelopes, preserving physical widths,
off-center Gaussians, callback conventions, interpolation and ramp reversal even
when evaluation duration changes. They acquire no inferred region metadata and
do not silently adopt new stretching rules. Old source calls must migrate
explicitly; loading old files is a separate compatibility path.

Native JLD2 callable persistence requires compatible definitions and trusted
files. Portable device bundles require registered callable codecs and explicit
opt-in for bundled Julia code. See [Persistence](@ref).

A runnable example is `docs/examples/pulses.jl`. Existing research notebooks are
not automatically rewritten during this breaking alpha API migration.
