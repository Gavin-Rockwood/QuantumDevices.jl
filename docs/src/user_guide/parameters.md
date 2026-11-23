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

Shared physical defaults can be declared with
`add_parameters!(device, Dict(EC1c=>4.3, EC2c=>4.5, EC12=>161.4))`.
Component and interaction declarations remain separate; `device.parameters`
returns a fresh assembled symbolic-keyed view. Use `add_parameters!` for updates
rather than mutating that view. Repeated equal defaults are permitted. Conflicting
defaults always fail with both values and sources, before recipe/build overrides
or controls are applied. Interaction defaults are optional until build, where
missing values can be supplied by explicit overrides. Recipe-local interaction
defaults stay local to that recipe.
