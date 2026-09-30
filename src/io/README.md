# Model bundles

```julia
save("my_device", model)   # creates a new directory; refuses existing paths
restored = load("my_device")
```

```
my_device/
  model.json                 # schema version, component definitions, interactions,
                             # truncations, dimension limit, gate index
  parameters.json            # all component and interaction parameter values
  components/                # optional custom-component JLD2 artifacts
    1.jld2
  gates/
    0001/
      gate.json              # duration, symbolic Hamiltonian, gate parameters
      pulses/                # optional generic-callable artifacts
        1.jld2
```

Gate directory numbers are storage identifiers. The model manifest preserves the
actual gate keys (including the distinction between strings and symbols), so keys
are never interpreted as filesystem paths. Gate parameters and built-in pulse
parameters are editable directly in `gate.json`. Component and interaction values
are editable in `parameters.json`; loading reconstructs the model's derived physics.
For tunable transmons, `EJ1` and `EJ2` are authoritative; `EJmax` and `d` are recomputed.
Custom components retain their stored operators and symbolic Hamiltonian.

Each structure owns its encoding and reconstruction rules in its own `.jl` file.
Pulse-specific rules live under `gates/pulses/`. Additional pulse types can extend
`_pulse_record(pulse, directory, index)` and
`_restore_pulse(::Val{:tag}, data, directory)`.

Built-in symbolic operations use an explicit allowlist. No stored source is parsed
or evaluated. Unsupported symbolic functions and parameter types fail explicitly.
Scalar integer, floating-point, Boolean, and complex values are supported; exact
Julia numeric types are not generally preserved by JSON.

Generic pulses store their callable in JLD2. Its defining module/type must be
available when loading in a fresh process. Arbitrary notebook closures are not
portable code archives; a reconstructed non-callable produces a load error.
Load bundles from trusted sources, particularly those containing JLD2 artifacts.

Saving writes a temporary sibling directory and moves it into place only after all
records succeed. Failed writes clean up that temporary directory. Existing bundles
are never overwritten. Derived model matrices, eigensystems, and tracking results
are rebuilt rather than persisted.

## Parameterized pulses

Both built-in and generic pulses expose `pulse.parameters` and evaluate `f(p, t)`.
`p.duration` is supplied by the gate; do not include it in stored parameters.

```julia
pulse = GenericPulseFunction((p, t) -> p.amplitude * sinpi(t / p.duration)^2,
                            (amplitude = 0.4,))
names = ["drive/amplitude"]
```

Calibration discovers finite real parameters identically for both pulse types.
Internal pulses store a registered name; generic pulses store their callable as a
JLD2 artifact and their parameters as JSON. Old `f(t)` / `f(t, duration)` generic
pulse records are rejected rather than silently reinterpreted.
