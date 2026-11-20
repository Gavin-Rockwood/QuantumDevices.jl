```@meta
EditURL = "tunable_coupler_control.jl"
```

# Tunable coupler control

Build two flux-tunable transmons and a coupler, then control the coupler's
actual Symbolics flux variable in the fixed reference basis. Frequencies use
GHz = E/h and time uses ns. `evolve` supplies the corresponding 2π factor.
This example constructs and evaluates Hamiltonians without running an optimizer.

```julia
using QuantumDevices, QuantumToolbox
import Symbolics

device = QuantumDevice(:tunable_coupler)
q1 = add_component!(device,make_flux_tunable_transmon(:q1,.198,11.143;dimension=31))
c = add_component!(device,make_flux_tunable_transmon(:coupler,.097,33.528;dimension=31))
q2 = add_component!(device,make_flux_tunable_transmon(:q2,.18,11.15;dimension=31))
recipe = Recipe(device,[q1,c,q2];states=Dict(q1=>0:3,c=>0:4,q2=>0:3))
```

Local aliases name actual Symbolics variables. Idle values live in separate
dictionaries; interacting components keep their existing variable identities.

```julia
EC1,ECc,EC2 = q1.variables["EC"],c.variables["EC"],q2.variables["EC"]
Symbolics.@variables EC1c EC2c EC12
add_parameters!(device, Dict(EC1c=>4.3,EC2c=>4.5,EC12=>161.4))
J1c = 8EC1*ECc/EC1c
J2c = 8EC2*ECc/EC2c
J12 = 8EC1*EC2*(1/EC12+ECc/(EC1c*EC2c))
n1,nc,n2 = recipe.inputs.q1.operators.charge,recipe.inputs.coupler.operators.charge,recipe.inputs.q2.operators.charge
for (name,coefficient,left,right) in (
    (:q1_coupler,J1c,n1,nc), (:q2_coupler,J2c,n2,nc), (:q1_q2,J12,n1,n2))
    add_interaction!(recipe,Interaction(name,coefficient*left*right))
end
model = build(device,recipe)
```

Runtime dictionaries use variables directly. Changing these coefficients does
not diagonalize a new reference basis; build another recipe for that purpose.

```julia
flux = c.variables["flux"]
H = make_hamiltonian(model,[flux];auto_time_dependence=false)
@assert H(Dict(flux=>0.),0.) ≈ make_hamiltonian(model)
H(Dict(flux=>.2),0.)
```

A control replaces the complete flux value. Choose its idle offset explicitly.

```julia
pulse = Pulse(GaussianZero(4.);duration=20.,amplitude=.1,offset=0.)
gate = build(device,GateRecipe(recipe;controls=[Control(:coupler_flux,flux,pulse)]))
@assert check_endpoints(gate;atol=1e-14)
make_hamiltonian(gate)(nothing,10.)
```

Propagation and calibration are explicit optional work. For example:
`result = evolve(gate;states=[(0,0,0),(1,0,0),(0,0,1),(1,0,1)])`.
Inspect `gate_matrix(result)` and leakage before accepting a
calibrated gate. See the gate recipe and calibration guides for those APIs.
