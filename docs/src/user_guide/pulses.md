# Pulses and flattop ramps

Pulses are finite scalar gate parameters evaluated as `pulse(t, duration)`.
`start` and `stop` use physical time, and `stop=nothing` uses the gate duration.
Built-in pulses return their `offset` outside the inclusive active window.

```@example pulses
using QuantumDevices, CairoMakie
constant = constant_pulse(0.2)
gaussian = gaussian_pulse(0.2, 0.15)
lobe = sine_squared_pulse(0.2, 0.5)
flat = ramped_flattop_pulse(0.2, 0.2)
@assert flat(0.0, 1.0) == 0
@assert flat(0.5, 1.0) == 0.2
fig = Figure(size=(700, 340))
ax = Axis(fig[1, 1], xlabel="Time", ylabel="Hamiltonian coefficient")
ts = range(-0.1, 1.1; length=301)
for (label, pulse) in (("constant", constant), ("Gaussian", gaussian),
                       ("sine squared", lobe), ("flattop", flat))
    lines!(ax, ts, [pulse(t, 1.0) for t in ts]; label)
end
axislegend(ax; position=:rt)
fig
```

A truncated Gaussian generally has nonzero values at its window endpoints.
It is not automatically shifted to a zero baseline. The periodic sine-squared
constructor takes frequency, not ramp duration; `frequency=0.5` gives one lobe
on normalized time `0:1`.

## One source for both ramps

```@example pulses
shared = ramped_flattop_pulse(0.2, 0.2;
    ramp=gaussian_pulse, ramp_kwargs=(; sigma=0.2))
@assert shared(0.0, 1.0) ≈ 0
@assert shared(1.0, 1.0) ≈ 0
shared(0.5, 1.0)
```

`ramp_kwargs` is a named tuple forwarded to the pulse constructor, with
`amplitude=1` unless overridden. You can also pass an existing pulse; kwargs then
override its named parameters. Shapes are resolved at construction, not at every
time sample. Named `:sine_squared`, `:linear`, and `:smoothstep` ramps remain
available through `available_ramps()` and do not accept kwargs.

## Independent rise and fall

```@example pulses
asymmetric = ramped_flattop_pulse(0.2, 0.2;
    ramp_up=gaussian_pulse,
    ramp_up_kwargs=(; sigma=0.2, center=0.4), split_up=0.4,
    ramp_down=sine_squared_pulse,
    ramp_down_kwargs=(; frequency=0.5), split_down=0.5)
@assert asymmetric(0.0, 1.0) ≈ 0
@assert asymmetric(0.2, 1.0) ≈ 0.2
@assert asymmetric(1.0, 1.0) ≈ 0
asymmetric(0.1, 1.0)
```

The shared `ramp`, `ramp_kwargs`, and `split` supply defaults for both sides;
`ramp_up`, `ramp_down`, and their respective kwargs/splits override them.
Rise and fall each take `ramp_time`; `2ramp_time` must fit in the active window.
Zero ramp time makes a constant window with abrupt edges.

## Source coordinates and normalization

Source pulses are evaluated with **duration 1** and time in `0:1`. Gaussian sigma,
center, start, and stop in ramp kwargs therefore refer to that normalized source,
not physical gate time. `split_up` selects the end of the rising source segment;
`split_down` selects the beginning of the falling segment. Both must lie strictly
between zero and one for pulse sources. Set splits to the intended peaks.

Each side subtracts its own endpoint value and divides by its endpoint-to-split
contrast. This guarantees the baseline and flat-top height for a suitable source,
but does not guarantee monotonicity or prevent overshoot if the source oscillates.
The split value must differ from the corresponding endpoint. Source amplitude and
offset cancel in this normalization; final height comes from flattop amplitude.

Pulse-source flattops are `GenericPulseFunction` objects. The resolved source
parameters live in the callable; `parameters(flattop)` exposes outer flattop
parameters, not a recursive list of ramp-source parameters. Reconstruct the pulse
to change source kwargs. See [pulse extension](../development/extensions.md) and
[API details](../reference/gates.md).
