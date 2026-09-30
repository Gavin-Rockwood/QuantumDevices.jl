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

## Parameterized controls

The simplest extension is `GenericPulseFunction(f, parameters)`, with `f(p,t)`.
`p.duration` comes from the owning gate and must not be stored in the parameters.

```@example extensions
pulse = GenericPulseFunction((p,t) -> p.amplitude * sinpi(t/p.duration)^2,
                             (; amplitude=0.2))
@assert pulse(0.5, 1.0) ≈ 0.2
parameters(pulse)["amplitude"]
```

This closure is suitable for interactive use. For persistence, define the callable
in a module as shown in the [fresh-process example](../user_guide/persistence.md).

## Custom pulse subtypes

```@example custompulse
using QuantumDevices
struct LevelPulse <: AbstractPulse
    value::Float64
end
QuantumDevices.pulse_value(pulse::LevelPulse, t, duration) = pulse.value
QuantumDevices.parameters(pulse::LevelPulse) =
    Dict{String,Tuple}("value" => ("value", pulse.value))
pulse = LevelPulse(0.2)
updated = setpath(pulse, "value", 0.4)
@assert pulse(0.0, 1.0) == 0.2
updated(0.0, 1.0)
```

The gate boundary verifies scalar values. A custom subtype may specialize
`QuantumDevices.validate_pulse(pulse, duration)` to check its domain and parameters;
that hook is qualified and is not exported. Parameter discovery must provide
paths that `getpath` resolves and `setpath` can reconstruct. The default structural
reconstruction uses the outer constructor; custom computed fields may need a
qualified `_replace_property` specialization.

## Calibration objectives

Pass `objective(model, candidate_gate, target)` to `calibration_problem`.
Define the full loss, subspace treatment, and leakage penalty explicitly. The
wrapper only discovers selected finite real parameters and reconstructs gates;
it does not choose optimizer algorithms, gradients, bounds policy, or noise metrics.

## Persistence hooks

Qualified internal hooks `QuantumDevices._pulse_record(pulse, directory, index)`
and `QuantumDevices._restore_pulse(::Val{:tag}, data, directory)` add pulse formats.
The record must contain a unique `"type"` tag; place artifacts inside the supplied
directory and use controlled relative filenames. These hooks are internal extension
points and may change with the bundle format. They are not a source-evaluation API.
Custom components use the existing component artifact path; supported symbolic
functions and parameter value types still constrain their records.

Here is a minimal JSON-only format for the `LevelPulse` defined above. In a real
extension, use a tag unique to your package and load these method definitions
before restoring the bundle.

```@example custompulse
QuantumDevices._pulse_record(pulse::LevelPulse, directory, index) =
    Dict("type" => "level_demo", "value" => pulse.value)
QuantumDevices._restore_pulse(::Val{:level_demo}, data, directory) =
    LevelPulse(Float64(data["value"]))
model = make_model([make_qubit("q", 0.0)], val(0), (;))
model.gates[:level] = DeviceGate((; drive=LevelPulse(0.2)),
                                param(:drive) * op(:q_x), 1.0)
mktempdir() do directory
    path = save(joinpath(directory, "device"), model)
    restored = load(path)
    @assert restored.gates[:level].parameters.drive isa LevelPulse
    @assert pulse_value(restored.gates[:level].parameters.drive, 0.5, 1.0) == 0.2
    println("Custom pulse format restored")
end
```
