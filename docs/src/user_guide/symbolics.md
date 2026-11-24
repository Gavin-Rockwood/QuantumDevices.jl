# Symbolic Hamiltonians

Scalar coefficients are ordinary Symbolics.jl variables and expressions. Quantum
operators are owner-aware `Sym` expressions, whose products preserve order.

```@example symbolics
using QuantumDevices
import Symbolics
import QuantumToolbox as QT
variables, idle = component_parameters(:q; EC=0.25, EJ=18.0)
EC, EJ = variables["EC"], variables["EJ"]
n, J = op(:q, :charge), op(:q, :jump)
expression = 4EC*n^2 - EJ*J/2
operators = Dict(OperatorRef(:q,:charge)=>QT.num(3)-1,
                 OperatorRef(:q,:jump)=>QT.destroy(3)+QT.destroy(3)')
defaults = Dict(variables[alias]=>value for (alias,value) in idle)
H = evaluate(expression, defaults, operators)
@assert H ≈ 4*0.25*operators[OperatorRef(:q,:charge)]^2 -
    18*operators[OperatorRef(:q,:jump)]/2
size(H)
```

Use scalar functions directly: `sqrt(EC)`, `cospi(flux)`, and products of scalar
variables require no wrapper. `val(number)` explicitly inserts a scalar leaf when
needed. Scalar arithmetic is delegated to Symbolics; it never turns operator
references into commuting scalar variables. In particular, `a*b` and `b*a`, and
`(a*b)'` with its reversed adjoint order, retain their quantum meanings. Scalar
functions of quantum operators are outside this sum/product compiler.

## Names, aliases, and idle values

`q.variables["EC"]` is an actual Symbolics variable. `q.parameters["EC"]`
is its numerical idle value. The string-keyed alias map is explicit; the implementation does
not infer local aliases from symbolic name strings. `component_parameters(name;
EC=.25, EJ=18.)` returns this pair of containers for custom components.

By default, variable names use the local alias followed by an underscore and
component name, for example `EC_q1` and `flux_coupler`. The helper and all builtin
component constructors accept `subscript`: `missing` uses the component name,
`nothing` produces names such as `EC`, and a nonempty Symbol or string supplies
an explicit suffix, for example `subscript=:left` produces `EC_left`. This is
independent of the component's display `label`.

Unsuffixed variables and shared explicit suffixes deliberately share symbolic
identities. Underscores can also make different alias/suffix pairs produce the
same name; choose different suffixes when sharing is not intended. Shared names
must have consistent idle defaults. Registration does not infer ownership from
name strings or reject deliberate sharing.
Identical component names and aliases recreate equal symbolic variables; distinct
components must have distinct names, and duplicate device component names are
rejected. Symbolics name equality, rather than object identity or metadata,
defines parameter identity. For intentionally shared external variables, construct
matching string-keyed variable and idle-value maps explicitly. Their names remain unchanged.
All idle values and overrides must be finite concrete numbers, and declaration
keys must be atomic scalar variables, rather than composite expressions.

## Interactions

```julia
g = Symbolics.variable(:shared_coupling)
interaction = Interaction(:exchange, g*a*b)
```

An interaction discovers every scalar variable used by its expression in
`interaction.variables`, keyed by full symbolic name. It stores the expression
and metadata. Numerical values belong to the consuming recipe or build:

```julia
add_interaction!(recipe,interaction)
recipe.parameters[g] = 0.02
model = build(device,recipe)
```

`device.variables` maps full symbolic names (for example, `"EC_q1"`) to
Symbolics variables. `device.parameters` groups each variable's resolved values
or `missing` by their consuming recipe and captured-history locations. Both are
fresh views; editing them does not update the device. Independent roots may use
different values. Conflicting assignments within one root report their values
and sources. Missing values are deferred until build and report required
variables and use sites. See the [parameter guide](parameters.md) for inventory
and assignment examples.

## Fixed-basis runtime values

`make_hamiltonian(model)` returns a static QuantumToolbox `QuantumObject`.
`make_hamiltonian(model, [flux])` returns a `QobjEvo` that reads `p[flux](t)`.
With `auto_time_dependence=false`, it reads `p[flux]` instead. Every other scalar
keeps its resolved build value. Runtime dictionaries use symbolic variables as
keys, including when multiple components have the same local alias. Numerical
conversion compiles scalar coefficients once, while keeping the model basis and
complete projected operator products fixed. Rebuild a recipe to recompute its
basis at another physical operating point.

`latex(expression; labels=Dict(...))` composes scalar mathematics with quantum
operator products. Symbolics owns its scalar display; this package does not
override global printing of Symbolics variables. Display labels never change
variable identity or supply idle values.

Recipe snapshots and captured calibration contexts preserve maps, aliases, and
symbolic name equality independently of later dictionary edits. Arbitrary callable
captured external state remains outside that boundary. Raw Symbolics structures
roundtrip with the installed JLD2/Symbolics versions; cross-version serialization
compatibility is not guaranteed. Complete device folders instead use explicit scalar AST records and registered operation IDs; unsupported metadata or operations require an explicit codec. Legacy standalone pulse persistence is unchanged.
