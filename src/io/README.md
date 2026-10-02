# Model bundles

```julia
save("my_device", model)   # creates a new directory; refuses existing paths
restored = load("my_device")
```

```
my_device/
  model.json                 # schema version, component definitions, interactions,
                             # writer's package release, truncations, dimension limit, gate index
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

`model.json` includes `quantumdevices_version`, taken from the loaded package's
version at save time, including prerelease suffixes. It describes the writer of
the most recent save, not an editing history. This is independent of
`schema_version`, which controls format compatibility. Model/component records
remain at `1`; gate records use `2` for timed recursive pulses. Legacy gates are
rejected clearly rather than silently reinterpreted. Missing
release metadata or a different writer release does not prevent loading a
supported schema. Future translators for breaking format changes can be placed
in `src/legacy/io/`; this change adds no migration infrastructure.

Each structure owns its encoding and reconstruction rules in its own `.jl` file.
Pulse-specific rules live under `gates/pulses/`. Additional pulse types can extend
`_control_record(control, directory, index)` and
`_restore_control(::Val{:tag}, data, directory)`.

Built-in symbolic operations use an explicit allowlist. No stored source is parsed
or evaluated. Unsupported symbolic functions and parameter types fail explicitly.
Scalar integer, floating-point, Boolean, and complex values are supported; exact
Julia numeric types are not generally preserved by JSON.

Custom `Envelope` and `Carrier` objects store their callable in JLD2.
Built-in pulses, envelopes, nested ramps, and carriers use JSON only. Its defining module/type must be
available when loading in a fresh process. Arbitrary notebook closures are not
portable code archives; a reconstructed non-callable produces a load error.
Load bundles from trusted sources, particularly those containing JLD2 artifacts.

Saving writes a temporary sibling directory and moves it into place only after all
records succeed. Failed writes clean up that temporary directory. Existing bundles
are never overwritten. Derived model matrices, eigensystems, and tracking results
are rebuilt rather than persisted.

## Timed pulses

```julia
pulse = Pulse(Gaussian(5.0); duration=20.0, amplitude=0.01,
    carrier=SineCarrier(4.6))
names = ["drive/amplitude", "drive/carrier/frequency", "drive/envelope/sigma"]
```

Pulses own duration and delay. Gate records preserve whether duration is explicit
or inferred; inferred timing is recomputed after loading and parameter updates.
Nested envelope, carrier, and custom named parameters remain selectable for
calibration. Old pulse APIs and records are intentionally unsupported.
