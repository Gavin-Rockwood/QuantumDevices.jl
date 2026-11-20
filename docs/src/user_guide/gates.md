# Gates and evolution

A pulse describes a coefficient; a `GateRecipe` says which physical variable it
controls and which extra Hamiltonian terms it drives. Build the recipe against a
device, evolve it, then project onto the states relevant to your gate.

## A complete X rotation

Start with a zero-idle qubit so that the control alone determines the motion.
With Hamiltonians in cycles per unit time, a constant `0.25 X` applied for one
time unit produces `exp(-i π X/2) = -i X`. The global phase does not affect gate loss.

```@example gates
using QuantumDevices, QuantumToolbox, CairoMakie
import Symbolics
device = QuantumDevice(:rotation)
q = add_component!(device,make_qubit(:q;frequency=0.0))
recipe = Recipe(device,[q])
A = Symbolics.variable(:drive_amplitude)
drive = Interaction(:drive,A*recipe.inputs.q.operators.x;parameters=Dict(A=>0.))
pulse = Pulse(Constant(;duration=1.);amplitude=.25)
definition = GateRecipe(recipe;terms=[drive],controls=[Control(:drive,A,pulse)])
prepared = build(device,definition)
result = evolve(prepared;saveat=range(0,1;length=41),abstol=1e-9,reltol=1e-9)
U = unitary(result)
target = ComplexF64[0 1;1 0]
@assert gate_infidelity(target,U) < 1e-7
gate_infidelity(target,U)
```

`evolve` supplies `2π` and merges pulse boundaries with solver stops. In a model
with nonzero idle Hamiltonian its drift remains present; there is no automatic
rotating frame or rotating-wave approximation. Saving on a dense grid requests
output samples; integration boundaries are supplied separately by `pulse_tstops`.

## Plot populations from the saved propagators

For a full propagator, choose input and output kets when projecting. Each column
is the evolved input; rows correspond to your chosen output basis.

```@example gates
basis_states = [basis(2,0),basis(2,1)]
matrices = [gate_matrix(result;states=basis_states,step=i).data for i in eachindex(result.times)]
ground = [abs2(M[1,1]) for M in matrices]
excited = [abs2(M[2,1]) for M in matrices]
fig = Figure(size=(700,320))
ax = Axis(fig[1,1];xlabel="Time",ylabel="Population",title="A constant-drive X rotation")
lines!(ax,result.times,ground;label="Initial state",color=:gray,linestyle=:dash)
lines!(ax,result.times,excited;label="Target state",color=:purple)
axislegend(ax; position=:lt)
ylims!(ax,0,1.05)
fig
```

The plot answers a population question for one input. A coherent gate also has
to act correctly on superpositions, which constrains relative phases. The default
`gate_infidelity(target,M)` includes leakage and those relative phases, while
ignoring a common global phase. See [conventions](../explanations/conventions.md).

## Evolve selected inputs

For a larger device, integrate only the columns needed for your target instead
of the full identity. This saves work without replacing the full Hamiltonian by
a smaller projected one.

```@example gates
columns = evolve(prepared;states=basis_states,abstol=1e-9,reltol=1e-9)
M = gate_matrix(columns)
@assert gate_infidelity(target,M) < 1e-7
leakage = 1 .- vec(sum(abs2,M.data;dims=1))
(M.data,leakage)
```

A selected-column result has no `unitary(result)`; it contains only its recorded
inputs. `gate_matrix` defaults to those same outputs, or accepts a distinct
`output_states` basis. Projected columns are never renormalized: lost norm is the
probability outside the output subspace. Ordered occupation labels can replace
kets for a model with tracked physical labels.

## Delays, offsets and endpoints

Controls replace a variable's full value, including its baseline. For a flux
excursion choose the idle offset explicitly. An additional charge-drive variable
can instead have a zero idle default supplied by its extra Interaction.

Gate duration defaults to the latest delayed control endpoint. An explicit
duration can include trailing idle time, but cannot cut a pulse short. Without
controls supply a duration. `check_endpoints(prepared;atol=...,rtol=...)` checks
matching with the captured reference explicitly; no smoothing or equality
constraint is imposed silently.

`make_hamiltonian(prepared)` returns the frequency-unit Hamiltonian. For direct
QuantumToolbox solver calls use `make_hamiltonian(prepared;scalar=2pi)` and include
pulse stops yourself. Built-in `evolve` already does both. `idle_phases(result)`
reports unwrapped reference idle phases without applying compensation.

Continue with [gate recipes](gate_recipes.md) for physical flux/microwave terms,
[calibration](calibration.md) for fitting controls, or the full
[transmon–resonator](../tutorials/transmon_resonator_control.md) and
[tunable-coupler](../tutorials/tunable_coupler_control.md) tutorials.
