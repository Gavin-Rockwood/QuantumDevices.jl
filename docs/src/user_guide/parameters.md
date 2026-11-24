# Parameter values and pulse updates

Physical scalar parameters are actual Symbolics.jl variables. A component exposes
local aliases through `q.variables`, while `q.parameters` stores numerical idle
values. For example, `q.variables["EC"]` is a symbolic variable and
`q.parameters["EC"]` is its idle charging energy. Composite scalar
expressions can appear in Hamiltonian coefficients. Component declaration maps
use local string keys; physical-default and override maps use atomic variables.
Values must be finite concrete numbers.

```julia
recipe = Recipe(device,[q]; parameters=Dict(q.variables["EC"]=>0.26))
model = build(device,recipe; parameters=Dict(q.variables["EC"]=>0.27))
@assert model.parameters[q.variables["EC"]] == 0.27
```

Build overrides take precedence over downstream recipe overrides, which take
precedence over upstream recipes and idle defaults. Shared variables must have
consistent declared idle values and resolved upstream values. Overrides never
alter the component's idle dictionary. The default component suffix keeps identical local aliases distinct. Explicit
shared suffixes or `subscript=nothing` can intentionally share variables.

`make_hamiltonian(model,[q.variables["EC"]])` reads runtime `p[q.variables["EC"]](t)`, or
`p[q.variables["EC"]]` with `auto_time_dependence=false`. Controls likewise target the
variable directly: `Control(:flux,q.variables["flux"],pulse)`. Control pulses replace
the complete value, including their explicitly chosen baseline offset.

## Editing saved recipes

Register a recipe with `add_recipe!` to capture its inputs, parameter overrides,
and local interactions. `draft(saved)` returns an independent editable tree:

```julia
trial = draft(saved)
trial.parameters[q.variables["EC"]] = 0.28
candidate = build(device,trial)
```

Snapshots and physical calibration contexts preserve aliases, parameter maps, and
symbolic name equality independently of later edits. Arbitrary callable captured
external state remains outside that snapshot boundary.

## Pulse paths

`parameters(pulse)` returns a cached flat `NamedTuple` of physical values, keyed by
symbols such as `Symbol("amplitude")`, `Symbol("carrier/frequency")`, and
`Symbol("envelope/sigma")`. Gate controls prefix paths with their control names,
for example `"drive/amplitude"`.

```@example parameters
using QuantumDevices
pulse = Pulse(Gaussian(; duration=1.0, sigma=0.1); amplitude=0.2)
updated = setpath(pulse,"amplitude",0.4)
@assert getpath(updated,"amplitude") == 0.4
@assert pulse.amplitude == 0.2
parameters(updated)
```

Paths use `/` separators and one-based indices. `setpath` copies the edited path
and rebuilds immutable structures; unchanged branches may remain shared. Use
`calibration_parameters` for selectable finite real optimization values, returned
as a flat path-to-`(path, value)` discovery map. Physical
Symbolics parameter maps are separate from this structural pulse-path interface.

## Inspecting device inputs

`device.variables` is a fresh dictionary from full symbolic names such as
`"EC_q1"` to atomic `Num` variables. `device.parameters` has the same keys;
each value is a dictionary grouping effective input values by their locations:

```julia
device.parameters["EC_q1"] # e.g. Dict(0.26 => [(source=:recipe, path=(:pair,), node=(:q1,))])
```

Locations identify the source category, repository path and recipe input path.
An unresolved use is grouped under `missing`; a declaration unused by any
inventoried root has an empty dictionary. Registered recipes, gates and supported
historical snapshots contribute usages independently. Captured historical values
are authoritative for those snapshots. Both views are independent read snapshots:
editing their dictionaries or location vectors does not change the device.

Component idle defaults remain on components. Interactions declare expressions,
symbols, and optional defaults. Recipe parameters and build overrides select
operating values. To
change a saved definition, edit a draft and register or replace it deliberately.
For a gate-only scalar, put its value in the gate's embedded model recipe draft,
leaving the registered base model unchanged. A missing static reference value
remains missing even if a control replaces that coefficient at runtime.


Interaction declarations may provide numerical defaults directly:

```julia
coupling = Interaction(:coupling, g * q1.operators["x"] * q2.operators["x"], Dict(g => 0.02))
add_interaction!(device, coupling)
device.defaults[string(g)] # 0.02
```

`device.defaults` combines registered component idle values and interaction defaults,
indexed by full symbolic names. Unassigned inputs appear as `missing`; conflicting
concrete declarations for the same symbol are rejected. This read snapshot is
separate from `device.parameters`, which groups recipe usages by resolved value.
Recipes look up each symbolic input by its full name in `device.defaults` and
capture the available interaction defaults when constructed, so later declaration edits
do not change existing recipes. Recipe overrides and explicit build parameters
take precedence over these captured defaults.
