# Saving and loading complete devices

`save_device` writes a shareable directory with a `device.yaml` manifest, a
content-deduplicated recipe DAG, gate records, and primitive numerical payloads.
`load_device` eagerly restores the registered device without building models,
calibrating gates, or evolving states. Captured calibration models are saved for
replay; an ordinary model or evolution is included only when explicitly attached.

```julia
add_gate!(device, (:CZ, :small_model), gate)
attach_artifact!(gate, :population_trace, evolution)
save_device("devices/my_device", device)
restored = load_device("devices/my_device")
gate = get_gate(restored, (:CZ, :small_model))
# Explicit work after loading:
# evolve(gate.calibration)
# problem = CalibrationProblem(gate.calibration)
```

The repository path is independent of `gate.name`. Duplicate insertion errors;
use `replace_gate!` explicitly. `move_gate!(device, :CZ, (:CZ, :small_model))`
changes the repository location without renaming the gate or changing its saved
calibration. Moves do no filesystem IO. Working recipes and historical snapshots
remain independent after loading even when their definitions share a saved recipe
node. Definition reuse does not cache numerical builds or merge tensor inputs.

Use `overwrite=true` to replace an existing directory. Saving preflights supported
records, writes a sibling staging directory, then renames it into place, retaining
the old destination for rollback until the replacement succeeds. Synchronize
concurrent edits yourself; saving is not a lock on mutable device data.

### Custom source and codecs

Unknown callables and metadata fail with their logical object path. Register a
stable codec ID/configuration rather than relying on anonymous closure capture:

```julia
module LabControls
using QuantumDevices
struct Width
    value::Float64
end
register_io_codec!("lab.width", Width; version=1,
    encode=x -> (; value=x.value), decode=c -> Width(c.value),
    source_id="lab_controls")
end
```

Supply its explicitly authored source with
`code_sources=Dict("lab_controls" => "LabControls.jl")`. For relative helper
includes use a descriptor with `root`, `entrypoint`, `module_name`, `files`, and
`dependencies`; dependencies are package names, optionally `Name=exact.version`. Missing packages or mismatched explicit versions fail without installation. Only declared contained files are bundled. Load with
`load_device(folder; allow_code=true)` to include the named module entrypoint in
`Main` and register its codecs automatically. This executes arbitrary Julia code
in your process; it is not a sandbox. No packages are installed automatically.
Identical registered module/source content is reused; conflicting content or an
unknown existing `Main` module requires a fresh session or another module name.
Loaded source bytes are retained so another save can re-bundle them after the
original folder moves.

### Portable results and optional native checkpoints

The default result stores an ordinary summary (`u`, native `retcode`, and objective
when available), while preserving owned `fitted_values`, all calibration settings,
parameterization, target/state order, validation and physical provenance. Replay
and continuation use those records and do not require the original optimizer.
The summary is not an `OptimizationSolution`; unsupported backend-specific data
is not silently converted into an objective or iteration count.

`save_device(...; native_checkpoint=true)` also saves original native optimizer
objects. Loading these requires both `load_native=true` and `trusted_native=true`
and matching Julia/package environment. Native deserialization is an independent
trust decision from bundled source execution. Portable loading never opens that
checkpoint. Legacy `save_pulse`/`load_pulse` remains available and retains its
previous native callable-serialization limitations; it is not the strict folder
format.

## Files and compatibility

The root `device.yaml` declares format/schema, registered gate paths and all
internal files with digests. `definitions/device.yaml` holds physical context;
`recipes/<content-id>.yaml` holds each distinct recipe node once;
`gates/<id>.yaml` holds working recipes and historical snapshots via references;
`payloads/arrays.jld2` contains raw primitive numerical bytes in little-endian order. Custom source trees
and optional native checkpoints occupy separate directories. IDs are local to a
save and independent of display names. Registered names may reference one recipe
node; embedded nodes do not pollute `device.recipes`.

Unknown schemas, codecs, operations, corrupt files, path traversal, symlinks and
unsupported metadata are rejected. Nothing is silently dropped or restored by
parsing printed symbolic expressions. Sparse/dense arrays, quantum dimensions,
construction-basis H/operators/products, tracked labels and selected propagation
columns are retained. Different resolved operating points keep independent
parameters/models even when their symbolic definition records are shared.

The schema is versioned separately from the package release. Existing native JLD
objects remain version-sensitive; old DeviceModel/DeviceGate folders are not
silently accepted or automatically migrated. `save_pulse(path,
result.recipe.controls[i].pulse)` saves only a pulse: it does not include the
calibration coordinates, link mapping, model provenance or duration policy.

FloquetBasis, FloquetSweep, AvoidedCrossingFit and ResonanceResult can also be attached explicitly. Their records contain data only; loading does not evolve a periodic drive, track modes or refit a crossing. A loaded basis requires an explicitly supplied Hamiltonian/params for new micromotion samples.

Pulse record version 2 distinguishes normalized Envelope/shape/region identifiers
from v1 physical-time shapes. V1 pulse records restore through bounded compatibility
adapters, preserving physical widths, off-center Gaussian positions and callback
conventions. They do not gain inferred regions or adopt normalized stretching.
Unknown versions and unexpected fields are rejected. New custom callbacks use
`f(parameters,u)` and portable bundles still require registered callable codecs.
