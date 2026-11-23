---
name: quantumdevices-controls
description: Run QuantumDevices pulse-to-evolution workflows and retrieve physical states, populations, energies, and gate matrices from prepared gates or saved results.
---

Start from existing artifacts: inspect a built model's state/energy directly; evolve an existing PreparedGate without rebuilding; retrieve from an existing GateEvolution without rerunning dynamics. For requested simulation, execute the necessary workflow within the requested scale. Inspection alone should not launch dynamics. Preserve notebook physical values, solver settings, run flags, outputs and unrelated edits.

## Complete tiny X rotation

Executable from the repository project; illustrative zero-idle qubit avoids drift. Hamiltonians use E/h in GHz, durations use ns. Constant `0.25 X` for 1 ns yields `-i X`.

```julia
using QuantumDevices, QuantumToolbox
import Symbolics

device = QuantumDevice(:rotation)
q = add_component!(device, make_qubit(:q; frequency=0.0))
recipe = Recipe(device, [q])
A = Symbolics.variable(:drive_amplitude)
drive = Interaction(:drive, A * recipe.inputs.q.operators.x;
    parameters=Dict(A => 0.0))
pulse = Pulse(Constant(; duration=1.0); amplitude=0.25)
gate_recipe = GateRecipe(recipe;
    terms=[drive], controls=[Control(:drive, A, pulse)])
prepared = build(device, gate_recipe)
model = prepared.model
H = make_hamiltonian(prepared) # QobjEvo in E/h; p is ignored by default
result = evolve(prepared; saveat=[0.0, 0.5, 1.0], abstol=1e-9, reltol=1e-9)
U = unitary(result) # final full propagator in construction coordinates
kets = [basis(2, 0), basis(2, 1)] # explicit coordinate basis for this X test
M = gate_matrix(result; states=kets)
ψfinal = U * kets[1]
p_excited = abs2(M.data[2, 1])
target = ComplexF64[0 1; 1 0]
@assert gate_infidelity(target, M) < 1e-7
@assert isapprox(p_excited, 1.0; atol=1e-7)
```

`GateRecipe` captures model definition, extra Interaction terms and named controls; `build` prepares its reference model and complete gate products. Extra terms do not enter the reference diagonalization. `result.times` includes the final duration; `step=i` selects a saved time in `unitary`, `gate_matrix` and `idle_phases`. Built-in `evolve` applies 2π and pulse boundary stops once; do not scale H separately for it. There is no automatic rotating frame, rotating-wave approximation or phase compensation.

## Retrieval and larger simulations

For an existing physical model use `ψ = state(model, (1,0))` and `E = stateenergy(model, (1,0))` with complete zero-based occupation tuples in `model.eigensystem.components` order. `model.energies[1]` is instead a one-based energy-order lookup. Nonzero-frequency qubit/transmon eigenstates need not equal coordinate `basis` kets.

For selected inputs, supply `evolve(prepared; states=[(0,0), (1,0)], saveat=...)` or ordered orthonormal kets. This propagates only those columns in the full model space. `gate_matrix(result)` then projects onto recorded inputs; `output_states=[...]` chooses another output basis. `unitary(result)` is unavailable for selected-column results. For a full result `gate_matrix(result)` defaults to energy-ordered reference eigenstates, so it can differ from `unitary(result)`.

For physical states use `unitary(result; step=i) * ψinitial` on full results; selected-column physical columns are in `result.propagators[i].data[:, j]`, ordered as supplied inputs, in construction coordinates. Populations in a chosen output basis are `abs2.(gate_matrix(result; ...).data[:, j])`. Projected columns retain leakage and are never renormalized. Reference `stateenergy` values and `idle_phases(result)` describe the captured idle model, not instantaneous driven energy expectation values. To request the latter, use `expect(H(nothing, t), ψt)` for a controlled QobjEvo and the evolved ket `ψt` in the same coordinates (or `expect(H, ψt)` for static H); the result retains E/h units.

## Controls and reuse contracts

- Controls replace the full target value: set the idle baseline with `offset=idle` (or atomic symbolic `offset=variable` resolved against the captured reference). They are not increments automatically added to idle.
- `Pulse(envelope, frequency; ...)` adds a cosine carrier; `Pulse(envelope; ...)` has none. Envelopes use physical local time. Delays shift pulse onset; sequence child delays are subsequent gaps. Duration defaults to the latest pulse end; without controls supply gate duration.
- `check_endpoints(prepared)` explicitly checks continuity with the reference; no smoothing is automatic. Choose pulse endpoints appropriate to the requested physics.
- Build once for repeated evolution of a prepared gate. `make_hamiltonian(model, gate_recipe)` compiles using a compatible existing reference without diagonalization; it requires all complete operator products already cached. A new uncached gate term needs an explicit `build(device, gate_recipe)`. Runtime coefficients retain the fixed projection basis.
- Direct QuantumToolbox solvers need `make_hamiltonian(prepared; scalar=2pi)` plus pulse stops. Built-in `evolve` already supplies both. Selected physical pulse paths can be supplied at runtime with `free_parameters`; see [make_hamiltonian](../quantumdevices-make-hamiltonian/SKILL.md).

Read only needed details: [pulses](../../docs/src/user_guide/pulses.md), [gate recipes/reuse](../../docs/src/user_guide/gate_recipes.md), [evolution/retrieval](../../docs/src/user_guide/gates.md). Use [Calibration](../quantumdevices-calibration/SKILL.md) for fitting physical controls and [Device I/O](../quantumdevices-device-io/SKILL.md) for requested persistence.
