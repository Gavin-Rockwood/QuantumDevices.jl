# Saving and loading model bundles

`save(path, model)` creates a **new directory** and returns its absolute path.
It refuses an existing destination. `load(path)` reconstructs physical derived
state and restores named gates.

```@example persistence
using QuantumDevices
model = make_model([make_qubit("q", 0.0)], val(0), (;))
model.gates[:X] = DeviceGate((; drive=Pulse(Constant(); amplitude=0.25, duration=1.0)),
                            param(:drive) * op(:q_x), 1.0)
mktempdir() do directory
    path = save(joinpath(directory, "device"), model)
    restored = load(path)
    @assert restored.H ≈ model.H
    @assert restored.gates[:X].duration == 1.0
    println(join(sort(readdir(path)), ", "))
end
```

## Bundle layout

```text
device/
  model.json          # schema, components, interaction, truncations, gate index
  parameters.json     # component and idle coupling values
  components/         # optional custom-component JLD2 artifacts
  gates/
    0001/
      gate.json       # duration, symbolic Hamiltonian, gate parameters
      pulses/         # optional generic-callable JLD2 artifacts
```

Numeric gate folders are storage identifiers, not gate names. The manifest retains
symbol versus string keys. Component and interaction values are editable in
`parameters.json`; gate values are in their gate records. Loading rebuilds
operators, the idle Hamiltonian, and dressed states. For tunable transmons, `EJ1`
and `EJ2` determine `EJmax` and `d` again. Custom components retain stored operators.

Saving stages a temporary sibling directory and moves it into place after all
records succeed. Failed saves clean up staging. There is no overwrite option.

## Release and format versions

`model.json` records `quantumdevices_version`, the version of the loaded
QuantumDevices package that performed the save, including any prerelease suffix
(for example, `0.1.0-alpha.2-DEV`). Each new save records its writer's version;
this field does not track the bundle's editing history.

`schema_version` separately identifies each record format. Model and component
records remain at `1`; gates use `2` for timed, recursive pulse records. Legacy
gate/pulse records are rejected with an explicit migration message.
Release metadata is informational: bundles without it remain readable, and a
different package release alone does not prevent loading a supported schema.
If a future format change breaks compatibility, schema translators can live in
`src/legacy/io/`. No migration infrastructure is introduced by this metadata field.

## Supported values and portability

Finite integer, floating-point, Boolean, complex, string, symbol, and `nothing`
values are supported by the JSON encoding. Exact Julia numeric types are not
generally preserved. Arbitrary nested parameter containers are not a supported
scalar value format. Symbolic operations use an explicit allowlist; unsupported
functions fail rather than being converted back into source code.

Built-in pulses, nested envelopes, ramps, and carriers are recursive JSON records
and reload in a fresh process without user definitions. Custom `Envelope` and
`Carrier` callables use JLD2 artifacts. Load their defining module/types before
loading a bundle. Notebook closures are not portable code archives. Use bundles
from trusted sources, especially those containing JLD2 artifacts.

## Saving a custom pulse shape

This step is only needed for custom pulse shapes. Built-in shapes such as
`Constant()` and `SineSquared()` already have their definitions in QuantumDevices.

A bundle stores a custom shape's data, but not its Julia code. Keep that code in
a file and include it before saving or loading the model.

**1. Define the shape in `PulseDefinitions.jl`.**

```julia
module PulseDefinitions
export SineLobe

struct SineLobe end
(::SineLobe)(p, t, duration) = sinpi(t / duration)^2
end
```

`SineLobe()` is an object that can be called like a function. `Envelope` calls it
with the pulse parameters `p`, time `t`, and pulse duration. This shape ignores
`p` and peaks at 1 halfway through the pulse.

```@raw html
<p>You can also download this file: <a href="../PulseDefinitions.jl" download>PulseDefinitions.jl</a>.</p>
```

**2. Create and save a model using that shape.**

Run this from the directory containing `PulseDefinitions.jl`. The destination
`device` must not already exist.

```julia
using QuantumDevices
include("PulseDefinitions.jl")
using .PulseDefinitions

model = make_model([make_qubit("q", 0.0)], val(0), (;))
pulse = Pulse(Envelope(SineLobe()); amplitude=0.2, duration=1.0)
model.gates[:shaped] = DeviceGate((; drive=pulse), param(:drive) * op(:q_x), 1.0)
save("device", model)
```

**3. Load it in a new Julia session.**

Use the same directory and load the shape definition before calling `load`.

```julia
using QuantumDevices
include("PulseDefinitions.jl")

restored = load("device")
restored.gates[:shaped].parameters.drive(0.5)  # approximately 0.2
```

Keep `PulseDefinitions.jl` with your application and share it alongside the
bundle. Saving a notebook closure alone does not make its code available in a
new session.