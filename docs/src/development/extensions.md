# Extending QuantumDevices

## Custom components

Provide matching string-keyed variable/default and operator/factory maps, and an owner-aware symbolic Hamiltonian. Constructors do not
generate operators or diagonalize a basis; the recipe builder does that work.

```@example extensions
using QuantumDevices
import QuantumToolbox as QT
variables, idle = component_parameters(:custom; gap=1.0)
operators, operator_generators = component_operators(:custom;
    z=dimension -> QT.sigmaz())
custom = Component(:custom; variables, parameters=idle, operators, operator_generators,
    dimension=2, hamiltonian=0.5variables["gap"] * operators["z"])
device = QuantumDevice(:example)
add_component!(device, custom)
model = build(device, Recipe(device,[custom]))
@assert length(model.spectrum)==2
custom.parameters["gap"]
```

To share a Symbolics variable across components, supply it unchanged in each
component's variable map and give it the same idle value under the corresponding local key in each parameter map.

## Custom envelopes

Use `Envelope(f, parameters; duration=T)` for a custom dimensionless
`f(p, u)` on normalized time `u∈[0,1]`. Amplitude, baseline, delay and carrier
belong to the enclosing `Pulse`; its duration forwards the envelope reference.

```@example extensions
shape = Envelope((p,u) -> p.scale * sinpi(u)^2, (; scale=1.0); duration=1.0)
pulse = Pulse(shape; amplitude=0.2, duration=1.0)
@assert pulse(0.5) ≈ 0.2
parameters(pulse)[Symbol("envelope/parameters/scale")]
```

Use the positional frequency overload `Pulse(shape, frequency; phase=0, reference=:pulse)`
for a constant carrier. Arbitrary custom carriers are not part of the current public API.
Interactive envelope closures work in memory; portable bundles need registered callable definitions.

## Custom envelope subtypes

```@example customenvelope
using QuantumDevices
import Symbolics
struct LevelEnvelope <: AbstractEnvelope
    value::Float64
end
QuantumDevices.envelope_value(shape::LevelEnvelope, t, duration) = shape.value
pulse = Pulse(LevelEnvelope(0.2); duration=1.0)
updated = setpath(pulse, "envelope/source/value", 0.4)
@assert pulse(0.0) == 0.2
updated(0.0)
```

Custom physical-time subtypes are adapted internally when a Pulse duration is supplied;
the example path descends through that adapter's `source` field. Prefer the common
normalized Envelope for new portable code. Custom subtypes implement
`envelope_value(shape, local_time, duration)`, and may
specialize `validate_envelope(shape, duration)`. Structural fields are discovered
recursively. Pulse evaluation requires finite
scalar values. The default immutable reconstruction calls the positional outer
constructor; computed fields may need a qualified `_replace_property` overload.

## Calibration objectives

Pass `objective(target_matrix, projected_operator)` to `CalibrationProblem`.
Define the full loss, subspace treatment, and leakage penalty explicitly. The
wrapper only discovers selected finite real parameters and reconstructs gates;
explicit algorithm arguments select native backends; the default is bound-aware NLopt Nelder–Mead. Gradients, bounds policy and noise metrics remain caller choices.

## Persistence hooks

Use `register_io_codec!` for explicitly registered custom types with versioned
encode/decode callbacks. See the [persistence guide](../user_guide/persistence.md)
for source registration and portability requirements. Builtin normalized envelopes
store shape configuration and reference duration; derived caches are reconstructed.
Do not rely on removed `_control_record` hooks or serialize anonymous closures as portable bundles.

Custom envelopes with joins or discontinuities can define
`envelope_tstops(shape, duration)` to return interior times relative to pulse
onset. `pulse_tstops` validates the stops, adds the pulse boundaries, and shifts
them by its delay. Carrier frequencies use cycles per unit time; IQ modulation
uses the convention documented in the [pulse guide](../user_guide/pulses.md).

## Tracking and labeling implementation

State tracking lives in `src/utils/tracking/`. Its `track_states.jl` loop matches
states, reorders their auxiliary data, updates continuity metrics, and retains
requested samples. Input checks live in `validation.jl`; matching methods live in
`trackers/quick.jl` and `trackers/global.jl`, each returning candidate indices and
matched amplitudes in reference order. `overloads.jl` handles array/property
behavior, while `display.jl` contains printing and key-range formatting.

All dressing code lives in `src/utils/tracking/dressed_states/`.
`dressed_states.jl` follows the product preparation, batched ramp diagonalization,
tracking, and final organization flow in one file, with options, result types,
and display methods alongside it. Model construction calls that feature.

Recipe validation has three boundaries. Constructors and edit methods report
immediate local mistakes. Capture then resolves names, checks active-path cycles,
and validates aliases, disjoint ownership, scoped operator bindings, current device
operator declarations, and selection forms while it copies the tree. Frozen recipe
storage does not certify compatibility with a different device. Device-less frozen
provenance may reuse storage; execution still passes through device-bound capture.

Preparation consumes that captured tree and resolves dimensions, index bounds,
limits, generator configuration, physical terms, and operating-point parameters.
Leaf generator configuration is retained in the disposable build node, so numerical
assembly takes no device argument and only checks actual generated QuantumObjects
and their dimensions. Occupation labels remain deferred until upstream models exist.
`make_hamiltonian` follows capture/preparation without allocating numerical operators
or enforcing the numerical tensor-size cap.

Default declarations are checked before operating-point overrides. Missing values
are diagnosed at the occurrence that needs them: an override on an unrelated input
branch cannot satisfy that use. Gate-only upstream assignments follow the owners
of their added operator terms and shared-variable values must agree across uses.
These checks do not introduce persistent validation flags, model caches, or additional
immutability guarantees.
