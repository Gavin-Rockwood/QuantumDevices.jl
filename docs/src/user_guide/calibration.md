# Gate calibration

`CalibrationProblem` captures a prepared reference model and ordered independent
coordinates: direct pulse paths, or explicitly declared Symbolics variables. `solve(problem, algorithm)` returns a `CalibrationResult`: its recipe
is ready to inspect or adopt, and its `.solution` is the untouched native optimizer
result. Wrapping performs no extra evolution or model build. `solve` is re-exported
from SciMLBase, so `using QuantumDevices` supplies the same function.

## Start with a small physical problem

This self-contained setup isolates the amplitude of the same constant-drive X
rotation used in the [evolution guide](gates.md). Construction is cheap; the solve
is shown as an explicit next action rather than run during documentation builds.

```@example calibration_setup
using QuantumDevices
import Symbolics
device = QuantumDevice(:fit_rotation)
q = add_component!(device,make_qubit(:q;frequency=0.))
recipe = Recipe(device,[q])
A = Symbolics.variable(:drive_amplitude)
term = Interaction(:drive,A*recipe.inputs.q.operators.x;parameters=Dict(A=>0.))
trial = GateRecipe(recipe;terms=[term],controls=[
    Control(:drive,A,Pulse(Constant(;duration=1.);amplitude=.20))])
prepared = build(device,trial)
target = ComplexF64[0 1;1 0]
problem = CalibrationProblem(prepared,["drive/amplitude"],target;
    states=[(0,),(1,)],lb=[.1],ub=[.4],
    evolution_settings=(;abstol=1e-9,reltol=1e-9))
problem.paths
```

```julia
result = solve(problem;maxiters=60) # default bound-aware NLopt Nelder–Mead
result.solution.retcode
checked = validate_calibration(result;save_evolution=true)
checked.validation.loss
accepted = QuantumGate(:X,checked)
```

The ideal amplitude is .25 in these units. Inspect the return code and an explicit
validation before adopting a result. The full tutorials extend this workflow to
[transmon rotations and a sideband](../tutorials/transmon_resonator_control.md)
and [coupler exchange](../tutorials/tunable_coupler_control.md), where leakage,
relative phases and basis convergence become central.

```julia
import SciMLBase
problem = CalibrationProblem(prepared, ["drive/amplitude"], target;
    states=[(0,), (1,)], lb=[0.0], ub=[0.5],
    evolution_settings=(; abstol=1e-9, reltol=1e-9))
result = solve(problem, algorithm; maxiters=100)
candidate = result.recipe
result.solution.retcode
result.solution.objective       # optimizer reported, not revalidated
```

The result retains ordered `coordinates`, owned `fitted_values`, affected physical
`paths`, concrete input `states`, the `target` matrix
in that same order, and `state_keys` when the inputs have complete occupation
labels. Explicit ket or mixed inputs have `state_keys=nothing`; their original
selection remains in `provenance.original_states`. Default labels are recorded
only when tracked reference rows prove their order. `virtual_z` is currently
`nothing`: no phase correction is fitted or applied.

## One calibration per working gate

```julia
gate = QuantumGate(:X, result)
gate.calibration === result
replay = evolve(gate.calibration)  # captured candidate and evolution settings
working = evolve(device, gate)    # current editable gate.recipe
```

`QuantumGate(name, recipe)` creates an uncalibrated working gate with
`calibration=nothing`. Adopting a result makes an independent working copy of its
recipe; editing that copy cannot change the result snapshot. A MaxIters/failed
solve with valid pulse values still produces a candidate result, but adoption
requires `QuantumGate(name, result; allow_unsuccessful=true)` explicitly. No
constructor runs a validation evolution or changes the native return code.

There is one singular calibration, with no named repository, implicit latest
selection, or named add/replace API. To replay a result, use `evolve(result)`;
there is no device-free `evolve(gate)` shortcut that might silently replay a stale
working recipe. Replay uses captured settings and inputs; explicit evolution
keywords override them. The current working recipe always uses the device path.

## Explicit validation

```julia
checked = validate_calibration(result;
    evolution_settings=(; abstol=1e-11, reltol=1e-11), save_evolution=true)
checked.validation.loss
checked.validation.evolution
@assert checked.solution === result.solution
@assert result.validation === nothing
```

Validation performs one explicit evolution, evaluates the captured objective and
target in the fixed state order, and returns a new independently captured result.
Its settings and optional saved evolution are under `.validation`. The original
result and optimizer-reported objective remain unchanged. Validation does not
silently change the fit's replay defaults or certify nonlinear constraints for
an arbitrary backend. Pulse constructors and calibration box feasibility remain
mandatory when constructing results.

For raw optimizer access, use `solve(problem.problem, algorithm; ...)`.
A native result can later be captured with `CalibrationResult(problem, raw)`;
`result.solution === raw`. Invalid/nonfinite/out-of-box optimizer parameters throw
on wrapping rather than fabricate a usable recipe. `calibrated_recipe(problem,
values)` remains available for external/raw parameter studies.

## Native optimization backends

`CalibrationProblem(...; adtype=SciMLBase.NoAD(), function_kwargs=(;))`
constructs a native `OptimizationFunction`. Supply native hooks such as `grad`,
`hess`, `hv`, `fg`, `cons`, and `cons_j` in `function_kwargs`; native problem
options such as `lcons/ucons` stay at construction level, while callbacks,
iteration limits and backend solve options belong to `solve`.

```julia
import OptimizationOptimJL
problem = CalibrationProblem(prepared, paths, target;
    function_kwargs=(; grad=gradient!), lb=lower, ub=upper,
    evolution_settings=(; abstol=1e-8, reltol=1e-8, remove_vacuum=true))
result = solve(problem,
    OptimizationOptimJL.Fminbox(OptimizationOptimJL.LBFGS()); maxiters=30)
```

`gradient!(g,x,p)` must be a valid derivative of this objective. Hessian and
constraint hooks keep their native SciML signatures. Load the chosen
Optimization-compatible solver package to activate its methods; there is no
QuantumDevices backend registration or hidden algorithm choice. Native results,
return codes and callbacks remain inside `.solution` unchanged.

Optional ADTypes choices pass through, for example
`ADTypes.AutoFiniteDiff(; absstep=1e-5, relstep=1e-5)` when the user's environment
loads ADTypes, FiniteDiff and its OptimizationBase extension. These are not new
mandatory dependencies. Configuration does not make the fixed complex
QobjEvo/ODE/recipe pipeline universally compatible with dual numbers or
reverse-mode AD. There is no silent differentiation fallback.

Raw objective probes may cross optimizer boxes without clipping when the pulse
is physically valid. `calibrated_recipe` and result construction still enforce
bounds; invalid durations, ramp geometry and nonfinite values still throw.
Choose finite-difference steps/directions within true physical domains and check
that inner ODE error is smaller than the differences being measured.

## Physical objective and energy origin

The default loss includes leakage:
`1 - (tr(M'M) + abs2(tr(T'M))) / (d*(d+1))`. `T` is unitary in the supplied
state order; `M` is the projected physical propagator, with columns never
renormalized. Relative/local/conditional phases count; global phase cancels.
A custom objective has signature `objective(target_matrix, projected_operator)`.

Use `evolution_settings=(; remove_vacuum=true, abstol=1e-4, reltol=1e-4)`
to retain a constant reference-vacuum energy subtraction in returned evolution.
It changes only the global origin and leaves default gate infidelity unchanged.
`energy_shift=:auto` is restored relative to that origin. Neither option is
forwarded to optimizer or native ODE-solver keywords. Target inputs belong in the
constructor's `states`, not inside `evolution_settings`.

## Continue a fit

```julia
next_problem = CalibrationProblem(result;
    evolution_settings=(; abstol=1e-8, reltol=1e-8))
next_result = solve(next_problem; maxiters=100)
```

This starts from the owned fitted coordinates and retained mapping, reusing the
captured physical model. It
retains the original state/target order, bounds and solve settings; evolution
overrides merge, preserving options such as `remove_vacuum`. Explicit
`function_kwargs` replace all derivative hooks. Validation settings do not become
warm-start defaults. Valid unsuccessful candidates can also seed another solve.

`CalibrationProblem(result, paths, target; ...)` explicitly selects direct physical
coordinates, even when the paths look unchanged. Explicitly supply `lb`, `ub`, and `function_kwargs` (use `nothing`, `nothing`,
and `(;)` to reset them). Integrality and constraint bounds reset unless supplied
again. Changing `adtype` with inherited hooks also requires explicit
`function_kwargs`. The caller must adapt callbacks when changing parameter or
target ordering; no mapping is guessed. Explicit `lb=nothing, ub=nothing` removes
inherited box bounds.

## Optimizer choice

`solve(problem; maxiters=60)` uses NLopt's native bound-aware Nelder–Mead
(`LN_NELDERMEAD`), without coordinate transforms. For this backend, `maxiters`
limits objective evaluations, rather than simplex iterations. It requires only
`using QuantumDevices`. Explicit `solve(problem, algorithm; ...)` selects any
compatible loaded SciML backend with that backend's own requirements and keyword
semantics. There is no public `calibrate` helper. Bounds and general constraints
are distinct; select an appropriate constrained backend when needed.

## Capture and persistence

`settings` captures evolution, derivative/function configuration and native
problem settings (including box/constraint bounds, integrality and sense).
`provenance.prepared` holds one reference model/recipe/component context for
replay; `provenance.objective` is the supplied metric. Avoid treating the cached
numerical model as the only persistent physical definition.

Arrays are ordinary owned data, not new read-only wrappers. Capture and gate
adoption copy user recipes, states, targets and mutable settings independently;
callables may still retain external mutable state. Native solutions can retain
problems, closures and backend-private objects. Native JLD round trips depend on
compatible code/dependencies; the complete-device folder format stores a portable summary plus owned settings, mappings and physical snapshots, with an optional trusted native checkpoint. No universal native solver-result migration is promised. Old QuantumGate/CalibrationRecord layouts and named storage APIs
must be rebuilt or manually migrated. Historical compensated losses are not
converted into physical fidelity.

## Symbolic pulse parameterization

Keep the gate's nominal pulses concrete. Declare the independent optimizer vector
and its physical outputs on CalibrationProblem:

```julia
Symbolics.@variables Aq Ac T
problem = CalibrationProblem(prepared, targetCZ;
    variables=[Aq => q1_excursion, Ac => coupler_excursion, T => q1_duration],
    links=Dict("q1_flux/amplitude" => Aq,
               "coupler_flux/amplitude" => Ac,
               "q1_flux/duration" => T,
               "coupler_flux/duration" => T-duration_gap),
    states=computational,
    lb=[0.11-q1_idle, 0.15-coupler_idle, 40.0],
    ub=[0.135-q1_idle, 0.25-coupler_idle, 120.0],
    evolution_settings=(; abstol=1e-6, reltol=1e-6, remove_vacuum=true))
result = solve(problem; maxiters=100)
refinement = CalibrationProblem(result;
    evolution_settings=(; abstol=1e-8, reltol=1e-8))
refined = solve(refinement; maxiters=200)
collect(zip(refined.coordinates, refined.fitted_values))
fitted_duration = getpath(refined.recipe, "q1_flux/duration")
```

Use an ordered vector or tuple of `variable => initial` pairs. Dictionaries are
rejected because native bounds, integrality, derivatives and `solution.u` use this
exact coordinate order. Variables are atomic Num identities compared by Symbolics
equality. The declaration's initial values materialize a captured starting candidate;
they may differ from nominal pulse values without modifying the supplied gate.
Expressions use only declared optimization variables and trusted scalar operations,
never physical model defaults, runtime time or other field paths. Registered functions
require their definitions/dependencies. Physical-symbolic offsets stay separate;
provide a concrete nominal offset if it should be optimized.

`problem.coordinates` has three entries in this example. `problem.paths` has four
canonical physical outputs, sorted for stable reporting in symbolic mode. They do
not define optimizer order or dimension. `parameterization` captures the immutable
mode/link map. Direct `CalibrationProblem(prepared, paths, target; ...)` remains
available, with `coordinates == paths` in caller order and no generated symbols.
`calibration_parameters(gate)` lists physical fields; there is no `parameters/T`
selection. Mixed direct-path and symbolic declarations, or half mappings, throw
before any keywords reach the native problem.

Bounds apply to coordinates. For gap 10.5 ns, delays 0/5 ns and two 5 ns ramps,
T in [40,120] gives a coupler duration T-10.5 with at least a 19.5 ns plateau.
Constructors reject invalid durations, ramps and interpolation nodes; native finite
differences can still encounter true geometry errors outside the box. No hidden
penalty or clipping is introduced. All mapped fields update atomically, including
nested sequences, interpolation nodes, carriers and custom Envelope parameters.
Handwritten gradients/Hessians/constraints need the full chain rule.

`result.fitted_values` is copied independently of mutable native `result.solution.u`.
It owns the continuation seed; `result.settings.problem.initial_values` records that
run's start. Chaining retains the mapping and coordinates without solving an inverse,
including nonlinear/noninvertible links. The concrete recipe and replay remain fixed
if native u changes later. Explicit new paths or variables/links request a coordinate
reset even when they look unchanged: supply lb, ub and function_kwargs explicitly;
integrality and constraint bounds reset unless supplied. Symbolic replacements must
supply variables and links together. Model replacement retains mappings/fitted values
and existing state-key/basis safeguards. Native solution identity, default/explicit
backends, validation and guarded QuantumGate adoption retain their existing contracts.
