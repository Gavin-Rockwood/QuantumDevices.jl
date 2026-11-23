# Voltage controls and static maps

A voltage-defined gate evaluates its pulse in volts, converts the voltage through
a captured static map to absolute normalized flux `Phi/Phi0`, and evaluates the
existing nonlinear Hamiltonian expression at that flux. Model reference flux
and waveform voltage baseline are independent choices.

```julia
q = make_flux_tunable_transmon(:q, 0.25, 12.0; dimension=15, flux=0.1)
device = QuantumDevice(:experiment)
add_component!(device, q)
map = AffineControlMap(:q_line; input=:voltage, output=q.variables["flux"],
    gain=0.4, offset=0.02)
add_control_map!(device, map)
recipe = Recipe(device, [q]; states=Dict(:q=>0:2))
voltage = Pulse(Gaussian(0.2; duration=1.0); amplitude=0.1, offset=0.2)
gate = GateRecipe(recipe; controls=[Control(:bias, ControlInput(:q_line, :voltage), voltage)])
prepared = build(device, gate)
H = make_hamiltonian(prepared; free_parameters=["bias/amplitude"])
H(Dict("bias/amplitude"=>0.15), 0.5)
```

`Recipe(...; control_maps=:auto)` captures eligible registry maps whose complete
output axes belong to selected components. Explicit `control_maps=[map]` or
`[:registered_name]` selects local additions or overrides. Upstream snapshots
retain their calibration; registry replacement affects subsequently constructed
recipes. Empty selection inherits upstream maps without registry additions.
`Recipe(device, upstream; ...)` creates a downstream recipe node. To attach a map
to the same physical endpoint, edit `draft(upstream).control_maps` as shown below.
`draft(recipe).control_maps` is editable. Captured maps appear at
`model.control_maps`; they never become physical defaults or optimization paths.
A gate's endpoint map is authoritative when reusing a physically compatible model.

Offsets in every mapped pulse, including nested sequences, must be explicit finite
real voltages. Known complex envelopes (DRAG, complex interpolation samples, or
complex segment scales) require a carrier, whose mixer produces a real voltage.
Custom envelope callbacks must return finite real values over the intended
waveform domain; this is a caller contract, not a continuous runtime check.
Sampled diagnostics can reveal violations at their supplied times.
Outside its active window the pulse's numeric voltage offset still
passes through the map. It need not reproduce the reference model flux.
`check_endpoints(prepared)` checks physical outputs against reference values.

For crosstalk use ordered `inputs`, `outputs`, a matrix `response`, and an output
`offset` vector. Rectangular and singular maps support forward evaluation. Each
selected map requires one control per input; hold an input with an ordinary
constant pulse. Direct physical controls and selected maps cannot overlap physical
outputs; two selected maps cannot overlap outputs either. Alternative maps may
coexist in a recipe until selected by controls.

```julia
map = AffineControlMap(:cross; inputs=(:left, :right),
    outputs=(q1.variables["flux"], q2.variables["flux"]),
    response=[0.4 0.03; -0.02 0.5], offset=[0.02, 0.01])
map_values(map, (0.1, 0.2)) # Ordered absolute flux tuple
```

Affine bounds are optional declared domains: `input_bounds` and `output_bounds`
contain ordered `(lower, upper)` pairs, or one pair for the scalar constructor.
`map_values` does not enforce affine bounds. `check_control_map(map, voltages)`
checks one point and throws on a failure. `control_map_diagnostics(gate, model;
times, params=nothing)` reports sampled input/output values and domain failures.
Finite samples cannot guarantee bounds over a continuous waveform. Values are
never clipped.

## Measured frequency to flux

Use the measured transition in the selected recipe model, including its coupled
components and approximation. The measurement callback explicitly selects the
transition and returns `E/h` frequency units, matching the experimental curve.
For time in ns, frequencies are GHz.

```julia
map = flux_map_from_frequency(device, response_recipe, frequency_of_voltage;
    name=:q_line, input=:voltage, parameter=q1.variables["flux"],
    measure=m -> m.energies((1,0,0)) - m.energies((0,0,0)),
    flux_interval=(0.0, 0.4), voltage_range=(-0.2, 0.8),
    frequency_tolerance=1e-4, measurement_label="q1 transition")
recipe = draft(response_recipe)
recipe.control_maps[map.name] = map
```

`frequency_of_voltage` can be a callable or ordered `(voltage, frequency)` samples.
Preparation builds the existing recipe at scanned physical flux values, holding
all other resolved physical references fixed. Choose an explicit sign/period and
strictly monotonic response branch. Constant or nonmonotonic physical response
and out-of-range measured frequencies error. A physical response span no larger
than `frequency_tolerance` is unresolved and also errors. Sweet-spot endpoints are permitted;
large `map.derivation.diagnostics.inverse_sensitivity` indicates poorly conditioned
flux inference. The measured voltage response itself can be nonmonotonic.
Applying one experimental curve to different model recipes produces distinct
effective maps for those approximations.

Preparation may perform many model builds. `flux_samples`, `voltage_samples`, and
`max_refinements` bound sampling/refinement; exhausting the budget before meeting
the requested sampled frequency tolerance errors. Held-out voltage midpoints use
actual model measurements, not just a comparison of interpolants. Diagnostics
record sampled residuals; they do not guarantee continuous or flux error bounds.
The internal explicit-response-table helper reports interpolation consistency
alone when no numerical response callback is supplied.

The returned scalar `TabulatedControlMap` uses piecewise-linear interpolation.
It stores pure knots and captured map-free physical preparation provenance in
`.derivation`; arbitrary experimental and measurement callables are not needed on
reload. A descriptive label does not recreate an arbitrary callback or promise
regeneration. No model build, inversion, or root finding occurs during evolution.

A manually prepared table uses the same forward path:

```julia
map = TabulatedControlMap(:q_line; input=:voltage, output=q.variables["flux"],
    voltage_knots=[-0.2, 0.0, 0.4], flux_knots=[0.02, 0.1, 0.3])
```

A table errors outside its voltage support (`extrapolation=:error`); choose support
covering waveform baselines and intended candidate ranges. No extrapolation or
clipping is available. Table knot crossings introduce slope changes. The solver
retains pulse timing stops; it does not root-find voltage-knot crossing times.
Check evolution tolerance convergence when these changes affect accuracy.

## Calibration and saving

Use existing pulse paths and links: `CalibrationProblem(prepared,
["bias/amplitude", "bias/delay", "bias/offset"], target; ...)`. Voltage parameters
remain ordinary scalar runtime parameters. Map coefficients stay fixed, and the
normal objective reuses the compiled Hamiltonian and pulse functions. Rebuild-only
envelope candidates also convert voltages before evaluating physical coefficients.
Accepted `result.recipe` retains its endpoint map snapshots.

`save_device` and `load_device` roundtrip affine/tabulated maps, measured derivation
data, mapped gates, and distinct registry/model/gate calibration snapshots. Folder
schema 2 reads schema 1 direct-only folders with empty map registries and empty
recipe maps. Native Julia struct checkpoints remain version-sensitive. Completed
builtin tables need no bundled executable calibration source. Custom static map
extensions need explicit compiler and IO support; a forward map does not imply an
inverse. Exact affine inverse waveform export is deferred.
