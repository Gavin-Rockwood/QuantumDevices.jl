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

## See the shape and its placement

Plot the envelope separately from its carrier when choosing widths and delays.
Here the same normalized Gaussian is stretched from 20 to 40 ns; the second
pulse also starts later. Neither change modifies the original envelope.

```@example pulse_shapes
using QuantumDevices, CairoMakie
e = Gaussian(;duration=20.,sigma=4.)
p = Pulse(e;amplitude=.2,delay=5.)
longer = Pulse(e;duration=40.,amplitude=.2,delay=10.)
t = range(0,55;length=401)
fig = Figure(size=(720,320))
ax = Axis(fig[1,1];xlabel="Time (ns)",ylabel="Control amplitude")
lines!(ax,t,p.(t);label="20 ns window, σ=4 ns")
lines!(ax,t,longer.(t);label="40 ns window, σ=8 ns")
axislegend(ax; position=:lt)
fig
```

A Gaussian has small nonzero endpoints. Choose `GaussianZero` when a baseline-
subtracted zero-ended window is intended; that is a different waveform. Delay
moves the complete pulse window, while amplitude scales its active envelope.

Compare a few analytic shapes at the same duration before choosing the window.
The Gaussian pair shows the effect of baseline subtraction; SineSquared and
Blackman taper smoothly, while RampedFlattop leaves a constant middle section.

```@example pulse_shapes
shapes = (
    ("Gaussian", Gaussian(;duration=20.,sigma=4.)),
    ("GaussianZero", GaussianZero(;duration=20.,sigma=4.)),
    ("SineSquared", SineSquared(;duration=20.)),
    ("Blackman", Blackman(;duration=20.)),
    ("RampedFlattop", RampedFlattop(4.;duration=20.)),
)
shape_times = range(0,20;length=401)
gallery = Figure(size=(900,280))
for (column, indices) in enumerate((1:2,3:4,5:5))
    ax = Axis(gallery[1,column];xlabel="Time (ns)",ylabel="Envelope")
    for index in indices
        label, envelope = shapes[index]
        lines!(ax,shape_times,envelope.(shape_times);label)
    end
    axislegend(ax; position=:lt)
end
gallery
```

## Regions and short compositions

`e.rise` / `e.fall` contain normalized interval metadata or `nothing`.
`e.rising` / `e.falling` select those actual regions, sharing the original payload.
The default view duration is the selected fraction of the reference duration.
Gaussian, SineSquared, Blackman and Slepian provide geometric half-window views;
these are not a guarantee of monotonicity. Missing region requests throw.
Public constructor intervals use physical times and are normalized once.

```@example pulse_sequence
using QuantumDevices, CairoMakie
window = SineSquared(; duration=10.)
e = EnvelopeSequence(
    (shape=window.rising,),
    (shape=Constant(; duration=20.),),
    (shape=window.falling,),
    (shape=Constant(; duration=2.), scale=0),
)
pulse = Pulse(e; amplitude=.2)
times = range(0,e.duration;length=401)
fig = Figure(size=(720,300))
ax = Axis(fig[1,1];xlabel="Time (ns)",ylabel="Control amplitude")
lines!(ax,times,pulse.(times))
vlines!(ax,[5.,25.,30.];color=:gray,linestyle=:dash)
fig
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

## Persistent recipes and compiled waveforms

`PulseRecipe(envelope; kwargs...)` owns the pulse definition; `Pulse(recipe)`
compiles it once. Existing `Pulse(envelope; kwargs...)` and
`PulseSequence(children...; kwargs...)` constructors remain available.
`PulseSequence` now returns a `Pulse` whose `recipe isa PulseSequenceRecipe`;
use recipe dispatch instead of testing `isa PulseSequence`.

`parameters(pulse)` returns its exact cached physical `NamedTuple`. Keys are flat
canonical names, including component envelope fields, for example
`parameters(pulse)[Symbol("envelope/segments/1/shape/sigma")]`.
`calibration_parameters(pulse)` retains the flat `path => (path,value)` discovery
contract. Segment starts and sequence durations follow component durations and
cannot be selected independently. Gate paths also accept segment envelope fields
without the `shape/` segment when unambiguous.

```julia
pulse = Pulse(Gaussian(2.; duration=10.); amplitude=0.3)
pulse(5.; amplitude=0.4)
pulse.f(5.; parameters(pulse)...)
pulse(5.; duration=20., envelope=(sigma=2.,))
```

Duration-only edits stretch physical widths, DRAG beta, and interpolation nodes.
An explicitly supplied width remains its supplied physical value, even when equal
to its original default. Flat parameter splats specify all supplied physical
fields explicitly. Unknown and structural keyword fields throw. `recipe`, `f`,
and cached `parameters` cannot be edited with `setpath`; physical edits return a
fresh compiled pulse.

Built-in sample data and explicit custom state are captured independently in
immutable storage. Custom callbacks retain `f(parameters,u)` with normalized
`u`; mutable opaque closure captures remain outside reproducibility guarantees.
Slepian bandwidth and custom-state changes require `setpath` and recompilation;
unchanged defaults can still be passed through `parameters(pulse)...`.
Standalone pulse records now store version-3 recipes and compile fresh callables
on load. Version-1 physical and version-2 normalized records remain readable.
