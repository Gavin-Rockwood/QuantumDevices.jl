# Extending QuantumDevices

## Custom components

Provide local quantum operators, named parameters, and a symbolic Hamiltonian.
The parent dimension must match every operator. A custom component's operators
are preserved during parameter updates. If changing its physical parameters also
changes its basis, register a reconstruction function in
`QuantumDevices.COMPONENT_CONSTRUCTORS` under the component's type string. The
constructor must accept `name`, `dimension`, and its stored parameters as keywords.
Registered components are reconstructed from parameters in saved bundles; unregistered custom
components use a JLD2 artifact instead. This registry is an internal extension
point and may change with the bundle format.

```@example extensions
using QuantumDevices, QuantumToolbox
custom = Component("custom", (; gap=1.0), (; z=sigmaz(), x=sigmax()),
                   0.5param(:gap) * op(:z), "custom", 2)
model = make_model([custom], val(0), (;))
@assert size(model.H) == (2,2)
model.parameters.custom_gap
```

## Custom envelopes and carriers

Use `Envelope(f, parameters)` for a custom dimensionless `f(p, t, duration)`
shape. Time is relative to pulse onset. Amplitude, baseline, and timing belong
to the enclosing `Pulse`.

```@example extensions
shape = Envelope((p,t,duration) -> p.scale * sinpi(t/duration)^2, (; scale=1.0))
pulse = Pulse(shape; amplitude=0.2, duration=1.0)
@assert pulse(0.5) ≈ 0.2
parameters(pulse)["envelope/scale"]
```

`Carrier(f, parameters; reference=:pulse)` accepts `f(p,t)` and chooses a pulse
or gate clock. Its named parameters are also discoverable. Interactive closures
work in memory; persistent custom callables need their defining code available
when loading. Built-in shapes and carriers require no user-defined module.

## Custom envelope subtypes

```@example customenvelope
using QuantumDevices
struct LevelEnvelope <: AbstractEnvelope
    value::Float64
end
QuantumDevices.envelope_value(shape::LevelEnvelope, t, duration) = shape.value
pulse = Pulse(LevelEnvelope(0.2); duration=1.0)
updated = setpath(pulse, "envelope/value", 0.4)
@assert pulse(0.0) == 0.2
updated(0.0)
```

Custom subtypes implement `envelope_value(shape, local_time, duration)`, and may
specialize `validate_envelope(shape, duration)`. Structural fields are discovered
recursively. Custom carrier subtypes implement
`carrier_value(carrier, local_time, gate_time)`. Pulse evaluation requires finite
scalar values. The default immutable reconstruction calls the positional outer
constructor; computed fields may need a qualified `_replace_property` overload.

## Calibration objectives

Pass `objective(model, candidate_gate, target)` to `CalibrationProblem`.
Define the full loss, subspace treatment, and leakage penalty explicitly. The
wrapper only discovers selected finite real parameters and reconstructs gates;
it does not choose optimizer algorithms, gradients, bounds policy, or noise metrics.

## Persistence hooks

Qualified internal hooks `QuantumDevices._control_record(control, directory, index)`
and `QuantumDevices._restore_control(::Val{:tag}, data, directory)` add envelope,
carrier, or pulse formats. Use a unique type tag and controlled local artifact
filenames. Load extension definitions before restoring their records.

```@example customenvelope
QuantumDevices._control_record(shape::LevelEnvelope, directory, index) =
    Dict("type" => "level_demo", "value" => shape.value)
QuantumDevices._restore_control(::Val{:level_demo}, data, directory) =
    LevelEnvelope(Float64(data["value"]))
model = make_model([make_qubit("q", 0.0)], val(0), (;))
model.gates[:level] = DeviceGate((; drive=pulse), param(:drive) * op(:q_x))
mktempdir() do directory
    path = save(joinpath(directory, "device"), model)
    restored = load(path)
    @assert restored.gates[:level].parameters.drive.envelope isa LevelEnvelope
    @assert restored.gates[:level].parameters.drive(0.5) == 0.2
    println("Custom envelope format restored")
end
```

These hooks are internal extension points and may change with the bundle format.
They never evaluate stored source. Component persistence and symbolic-operation
allowlists continue to apply independently.
