# Saving and loading model bundles

`save(path, model)` creates a **new directory** and returns its absolute path.
It refuses an existing destination. `load(path)` reconstructs physical derived
state and restores named gates.

```@example persistence
using QuantumDevices
model = make_model([make_qubit("q", 0.0)], val(0), (;))
model.gates[:X] = DeviceGate((; drive=constant_pulse(π/2)),
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

## Supported values and portability

Finite integer, floating-point, Boolean, complex, string, symbol, and `nothing`
values are supported by the JSON encoding. Exact Julia numeric types are not
generally preserved. Arbitrary nested parameter containers are not a supported
scalar value format. Symbolic operations use an explicit allowlist; unsupported
functions fail rather than being converted back into source code.

Internal pulses store a name and parameters. Generic pulses and pulse-source
flattops store callable artifacts in JLD2. Load their defining module/types before
loading a bundle. Notebook closures are not portable code archives. Use bundles
from trusted sources, especially those containing JLD2 artifacts.

## A callable restored in a fresh process

```@raw html
<p>Download the module-owned callable definition: <a href="../PulseDefinitions.jl" download>PulseDefinitions.jl</a>.</p>
```

The module defines a callable type. In an application, keep this module in your own codebase
and include or import it before loading the bundle.

```@example portable
using QuantumDevices
definitions = joinpath(pkgdir(QuantumDevices), "docs", "src", "public", "PulseDefinitions.jl")
Base.include(Main, definitions)
model = make_model([make_qubit("q", 0.0)], val(0), (;))
pulse = GenericPulseFunction(Main.DocsPulseDefinitions.SineLobe(), (; amplitude=0.2))
model.gates[:shaped] = DeviceGate((; drive=pulse), param(:drive) * op(:q_x), 1.0)
mktempdir() do directory
    path = save(joinpath(directory, "device"), model)
    environment = dirname(Base.active_project())
    script = """
    using QuantumDevices
    include($(repr(definitions)))
    m = load($(repr(path)))
    @assert pulse_value(m.gates[:shaped].parameters.drive, 0.5, 1.0) ≈ 0.2
    println("Fresh-process restore passed")
    """
    print(read(`$(Base.julia_cmd()) --startup-file=no --project=$environment -e $script`, String))
end
```
