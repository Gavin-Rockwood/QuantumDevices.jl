# # Transmon–resonator control
#
# Scalar Hamiltonian coefficients use Symbolics expressions; operator products
# preserve their order and are projected as complete products. Frequencies are
# GHz = E/h and times are ns. This example performs no optimization.

using QuantumDevices, QuantumToolbox
import Symbolics

device = QuantumDevice(:transmon_resonator)
t = add_component!(device,make_transmon(:transmon,.25,18.;dimension=31))
r = add_component!(device,make_resonator(:resonator;frequency=6.0,dimension=5))
recipe = Recipe(device,[t,r];states=Dict(t=>0:3,r=>0:3))
n = recipe.inputs.transmon.operators.charge
a = recipe.inputs.resonator.operators.a
g = Symbolics.variable(:charge_coupling)
add_interaction!(recipe,Interaction(:coupling,g*n*(a+a'),Dict(g=>.03)))
model = build(device,recipe)

# Rebuild with an explicit variable-keyed override to change the physical basis.
changed = build(device,recipe;parameters=Dict(t.variables["EJ"]=>19.))
@assert changed.parameters[t.variables["EJ"]] == 19.
@assert t.parameters["EJ"] == 18.

# A microwave amplitude need not occur in the idle Hamiltonian. Declare it on
# an extra gate term; its operator is projected in the existing reference basis.
amplitude = Symbolics.variable(:microwave_amplitude)
drive = Interaction(:charge_drive,amplitude*n,Dict(amplitude=>0.))
pulse = Pulse(GaussianZero(4.),5.0;duration=20.,amplitude=.02)
gate = build(device,GateRecipe(recipe;terms=[drive],
    controls=[Control(:drive,amplitude,pulse)]))
H = make_hamiltonian(gate)
H(nothing,10.)

# Controls use captured pulses rather than runtime dictionaries. To alter a pulse,
# create a draft or reconstruct a GateRecipe and build it again. Physical evolution
# is optional: `evolve(gate;states=[(0,0),(1,0)])` applies 2π once and preserves
# projected column norms, so leakage remains measurable.

# ## Optional periodic sideband and Stark diagnostic
#
# This scans a cosine plateau, not the finite Gaussian gate above. Physical
# component operators avoid the original recipe's now-captured input scope.
# Inspect fit residuals, tracking overlaps and competing-mode flags; a resonance
# estimate does not validate a shaped pulse or its two-state approximation.
run_resonance = false
if run_resonance
    sideband_states = [(2,0), (0,1)]
    bare = abs(stateenergy(model,sideband_states[2])-stateenergy(model,sideband_states[1]))
    frequencies = range(0.9bare,1.1bare;length=17)
    resonance = find_resonance(model,t.operators["charge"],0.02,frequencies;
        states=sideband_states,phase=0,
        evolution_settings=(;abstol=1e-10,reltol=1e-10))
    shift = stark_shift(resonance)
    figure, axes = plot_resonance(resonance;frequency_offset=bare)
end
# To retain a former sine drive, set phase=-pi/2 explicitly.
