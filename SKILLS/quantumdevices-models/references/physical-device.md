# Two tunable transmons and a tunable coupler

Complete construction example, using illustrative parameters and deliberately small charge cutoffs for an inexpensive starting point. Preserve the user's device values and determine convergence before quantitative use. `EC`, `EJ`/`EJmax` and coupling coefficients use E/h in GHz; flux is Φext/Φ0, asymmetry `d` and offset charge `ng` are dimensionless. Charge-basis dimensions must be positive odd counts (2*ncut+1).

```julia
using QuantumDevices

device = QuantumDevice(:two_qubits_coupler)
q1 = add_component!(device, make_flux_tunable_transmon(:q1, 0.27, 17.0; d=0.1, flux=0.0, dimension=9))
c = add_component!(device, make_flux_tunable_transmon(:c, 0.30, 25.0;
    d=0.1, flux=0.2, dimension=9))
q2 = add_component!(device, make_flux_tunable_transmon(:q2, 0.29, 19.0; d=0.1, flux=0.0, dimension=9))
for (name, a, b) in ((:q1c, q1, c), (:cq2, c, q2))
    add_interaction!(device, Interaction(name,
        0.02 * a.operators["charge"] * b.operators["charge"]))
end
recipe = Recipe(device, [q1, c, q2];
    states=Dict(q1 => 0:2, c => 0:2, q2 => 0:2), interactions=:auto)
flux = c.variables["flux"]
model = build(device, recipe; parameters=Dict(flux => 0.23))
H = make_hamiltonian(model)
ψ100 = state(model, (1, 0, 0))
E100 = stateenergy(model, (1, 0, 0))
@assert size(H) == (27, 27)
@assert model.eigensystem.components == (:q1, :c, :q2)
```

For a flux excursion, define a numerical pulse with its full baseline, e.g. `Pulse(RampedFlattop(1.0; duration=8.0); amplitude=0.01, offset=0.23)`, and bind `Control(:coupler_flux, flux, pulse)` in `GateRecipe(recipe; controls=[...])`. Preserve the build override when preparing: `prepared = build(device, gate_recipe; parameters=Dict(flux => 0.23))`; an unqualified rebuild from `recipe` would return to its 0.2 default. For compilation alone, `make_hamiltonian(model, gate_recipe)` reuses this captured 0.23 reference and cached flux products without rebuilding. Then follow [Controls](../../quantumdevices-controls/SKILL.md). An atomic symbolic `offset=flux` instead captures the prepared model's reference flux during compilation. These are illustrative controls, not an implemented calibrated entangling gate.

The built-in effective Josephson energy is `EJmax * sqrt(cospi(flux)^2 + d^2*sinpi(flux)^2)`; do not substitute an independent EJ parameter for this component. Verify constructor details in [source](../../../src/Components/constructors/flux_tunable_transmon.jl). For larger grouped reductions read the [quickstart](../../../docs/src/getting_started/quickstart.md); do not duplicate upstream internal couplings.

For fixed-frequency qubits replace q1/q2 constructors with `make_transmon(name, EC, EJ; dimension=9)`; preserve the intended physical model.
