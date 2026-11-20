# Gate recipes

`GateRecipe` binds named controls to an ordinary model recipe and optionally adds
symbolic Hamiltonian terms. The executable [three-transmon example](https://github.com/Gavin-Rockwood/QuantumDevices.jl/blob/dev/docs/examples/gate_recipes.jl)
uses the two-qubit and tunable-coupler setup from `dev/devV7.ipynb`: first a flux
excursion, then a microwave charge drive, then both together. Run it from the
repository root with `julia --project docs/examples/gate_recipes.jl`.

```julia
flux_control = Control(:coupler_flux, ϕc, flux_pulse)
Ω = Symbolics.variable(:microwave_amplitude)
drive_term = Interaction(:q1_drive, (Ω) * recipe.inputs.q1.operators.charge;
    parameters=Dict(Ω=>0.0))
recipe_with_controls = GateRecipe(recipe;
    controls=[flux_control, Control(:q1_microwave, Ω, microwave_pulse)],
    terms=[drive_term], duration=24)
prepared = build(device, recipe_with_controls)
H = make_hamiltonian(prepared)
H(nothing, 10.0)
```

Each control replaces its target's full value. Set the pulse or sequence offset
to the desired baseline; it is not automatically added to the idle value. Named
controls must have distinct names and targets. Compose pulse contributions inside
one pulse or envelope when they share a target. Uncontrolled parameters keep
resolved recipe values, with explicit build parameters taking precedence.
Runtime `p` in the returned `QobjEvo` is ignored: pulses and fixed values are
captured at build time. Without controls, `make_hamiltonian(prepared)` returns a static
`QuantumObject`.

Duration defaults to the latest control end, including delays and sequence gaps.
An explicit duration can add trailing idle time but cannot cut off a control.
Without controls, supply a duration. No endpoint smoothing or equality constraint
is automatic. `check_endpoints(prepared; atol=..., rtol=...)` explicitly checks
that every control matches its reference value at time zero and the gate end.
`pulse_tstops(prepared)` supplies the union of pulse boundaries and gate endpoints
for future solver use.

Extra terms use `Interaction` and discover every coefficient variable, including
variables absent from the physical model Hamiltonian. Defaults may be supplied
by the term, device, or explicit recipe/build operating-point assignments.
Conflicting defaults are rejected before overrides or controls are applied. Scoped
operators must belong to the endpoint recipe's own exposed inputs. Physical
operator references also work, including on a `Leaf` endpoint. Products are
projected as complete products before truncation: `P' * A * B * P` is preserved.
Extra terms never enter the reference diagonalization; `prepared.model` retains
the ordinary model spectrum, basis, and state labels. `prepared.terms` contains
the physical symbolic extra terms; their projected matrices live in the endpoint
model's operator-product cache. Model and gate preparation share default-source
validation, operating-point precedence, and complete-product collection. Each
required local product is projected once before full words are tensored; reference
H, interactions, exposed operators, and compilation reuse that dictionary. No
upstream numerical models are retained.

The build captures an independent recipe tree, controls, and term definitions.
Use `draft(prepared.recipe)` to obtain an editable copy, or reconstruct a new
`GateRecipe` with changed controls. This is not a calibration snapshot of device
component definitions. Arbitrary callable captured state also remains outside
the snapshot boundary.

`make_hamiltonian` defaults to `scalar=1`, preserving `E/h` frequency units. The static
model, free-parameter model, prepared gate, and cached model plus `GateRecipe` or
`QuantumGate` forms all accept `scalar::Number`. Use `scalar=2pi` for direct
QuantumToolbox evolution in GHz and ns. The scalar and fixed multiplicative
coefficient factors are folded once into fresh matrices; nonlinear and additive
expressions remain in runtime coefficients. `herm_check` applies to the final
scaled matrices. The cached model operators are unchanged. Built-in `evolve`
applies its own `2π` conversion, so do not add a second factor.

See the [gate API reference](../reference/gates.md) for constructor signatures.

## Physical gate evolution

`evolve(prepared)` propagates the full identity and saves only the final physical
lab-frame unitary in the model construction coordinates. `unitary(result)` retrieves it.
`model.states(key)` and `state(model, key)` are valid eigenstate inputs in those
same coordinates; no separate diagonal-basis conversion is needed. `evolve` applies the `2π`
frequency conversion once, resolves pulse boundaries, and forwards native
QuantumToolbox solver keywords such as `abstol`, `reltol`, and `dtmax`.

For selected computational inputs, use ordered occupation labels in the reference
model's component order (or orthonormal QuantumObject kets):

```julia
result = evolve(prepared; states=[(0,0,0), (1,0,0)], saveat=[0, 2, 4])
physical = gate_matrix(result)
phases = idle_phases(result) # labeled, unwrapped 2π E t
```

The final duration is always saved, including when omitted from `saveat`.
`step=3` selects a saved sample for all three retrieval functions. Selected-input
results store physical evolved columns and cannot supply a full unitary.
Their input order is fixed; `output_states` may specify a distinct ordered
orthonormal projection basis. Full results accept `states` for their input basis. With no selection,
`gate_matrix` uses energy-ordered reference eigenstates for inputs and outputs;
this logical gate matrix can differ from the construction-coordinate `unitary`.

All local and conditional dynamical phases remain physical. Automatic energy
centering subtracts a constant only during integration and restores its global
phase exactly. Projected columns retain leakage and are never renormalized.
Reference energies, labels, component ordering, and provenance are
captured independently of later model edits. Arbitrary mutable state hidden
inside user callbacks or pulse closures is outside this snapshot guarantee.
No corrected propagators or phase histories are stored. Solver failures throw.

## Reusing a reference model

`make_hamiltonian(model, gate)` compiles a `GateRecipe` (or a `QuantumGate`'s current
working recipe) using an existing model's cached complete operator products:

```julia
flux = component.variables["flux"]
recipe = Recipe(device, [component]; states=Dict(component => 0:3))
model = build(device, recipe; parameters=Dict(flux => 0.2))
pulse = Pulse(GaussianZero(1; duration=4); duration=4, amplitude=0.03, offset=flux)
gate = GateRecipe(recipe; controls=[Control(:flux, flux, pulse)])
hamiltonian = make_hamiltonian(model, gate)
```

The atomic Symbolics `offset=flux` reads the captured `model.parameters[flux]`
once during compilation. It stays symbolic in the gate definition and does not
read the runtime controlled value. Numeric offsets retain their ordinary behavior.
The same gate can therefore use another model's reference operating point when
its recipe leaves that assignment unspecified. Sequence offsets and child offsets
follow the same rule, including gaps and times outside the pulse. An unresolved
pulse cannot be called directly; compile it against a model first. Symbolic offsets
appear in parameter discovery but cannot be selected as numeric optimizer variables.

Reuse compares the concrete recipe tree, component/input order, state selections,
generator choices, interactions, and dressing choices with captured provenance.
Independent equivalent scopes and snapshots are supported. Effective explicit
parameter assignments must match the model, with downstream assignments overriding
upstream assignments; omitted assignments inherit captured values. Execution thread
settings do not affect compatibility. Named roots or upstream references require a
concrete captured tree; this overload never consults a live device registry.

Compilation performs no diagonalization, state tracking, generator calls, or
implicit build. Extra terms must use products already present in
`model.operator_products`; missing products require an explicit `build(device, gate)`.
Products preserve their full order and adjoints through projection. Extra defaults
must agree with existing model reference values before controls replace coefficients.
`check_endpoints(model, gate)` uses these same resolved pulse references. Hamiltonians
retain `E/h` units without a `2π` factor, returning a `QobjEvo` with controls and a
`QuantumObject` without controls.


Added terms share the model's symbolic parameter identities. If a term uses a
parameter already resolved in the reference Hamiltonian, that captured value is
authoritative: applicable explicit recipe assignments must agree with it. For
example, an otherwise-unused assignment on the target operator's input cannot
change another component's reference frequency during gate compilation. Without
a conflicting assignment, using that same parameter on another owner's operator
legitimately uses the shared reference value. Controls remain the explicit mechanism
for replacing a reconciled parameter value at runtime.

Genuinely additional gate variables resolve from applicable owner-scoped recipe
assignments, then idle defaults; all terms using the same identity must agree.
Fresh builds honor explicit build/downstream overrides before this reconciliation.
Cached conversion checks its gate recipe against captured reference values and
requires compatible declared operating points. Conflicting default declarations
are still rejected before operating-point overrides are applied.

`gate_matrix` always projects the physical propagator, preserving local and
conditional dynamical phases. `idle_phases` is a reference-energy diagnostic;
it does not change the matrix. A global integration energy shift is restored
exactly before results are returned.

`evolve(prepared; remove_vacuum=true)` retains a constant reference-vacuum energy
subtraction in its returned propagators. It changes only the global energy origin;
`energy_shift` still restores its own integration-centering phase. The result's
solver metadata records `remove_vacuum`, `vacuum_energy`, and `energy_shift`.
Neither energy-origin option is forwarded to the native ODE solver.

## Numerical pulse geometry and optimization

GateRecipe owns physical pulses, extra terms, metadata and explicit duration policy.
Optimization variables and symbolic pulse links belong to CalibrationProblem; see
[calibration parameterization](calibration.md#Symbolic-pulse-parameterization).
Ordinary `Pulse` and custom `Envelope(f, parameters)` constructors receive numerical
geometry. Existing atomic Symbolics offsets still refer to physical model defaults.

An omitted gate duration follows the latest endpoint after a numerical update.
An explicitly supplied duration stays fixed even when it originally equals that
endpoint. Sequence durations derive from segment/child geometry. Physical-time
interpolation nodes must be updated alongside duration; no automatic stretching or
clipping occurs. The calibration updater applies all mapped fields together before
constructor validation, including custom Envelope parameters and carrier fields.

Use `getpath(gate, "drive/duration")` or existing structural paths such as
`controls/1/pulse/duration` to inspect the concrete candidate. `setpath` edits a
physical field and returns a new physical gate; it does not change a calibration
mapping captured on a result. There is no optimization parameter namespace or
`setparameters` API on GateRecipe.

## Inspecting gates and calibration

`display(gate)`, `display(problem)` and `display(result)` show stored summaries:
physical control blocks, optimizer coordinates/bounds, and captured fit/validation
status. `print` uses a bounded one-line summary. A QuantumGate shows its current
working recipe separately from its saved calibration. Printing performs no build,
evolution or validation; `result.solution` retains the backend's own display.
