# Gate recipes: two flux-tunable qubits and a tunable coupler.
# Run with `julia --project docs/examples/gate_recipes.jl`.
# Frequencies are GHz = E/h and times are ns. Hamiltonian evaluation below
# introduces no 2π; evolve applies that factor once during integration.
using QuantumDevices
import Symbolics
using QuantumToolbox

device = QuantumDevice(:two_qubits_and_coupler)
q1 = add_component!(device, make_flux_tunable_transmon(:q1, 0.25, 18.0; d=0.1))
coupler = add_component!(device, make_flux_tunable_transmon(:coupler, 0.30, 30.0; d=0.1))
q2 = add_component!(device, make_flux_tunable_transmon(:q2, 0.26, 19.0; d=0.1))
ϕ1, ϕc, ϕ2 = q1.variables["flux"], coupler.variables["flux"], q2.variables["flux"]
EC1, ECc, EC2 = q1.variables["EC"], coupler.variables["EC"], q2.variables["EC"]
Symbolics.@variables EC1c EC2c EC12
recipe = Recipe(device, [q1,coupler,q2];
    initial_dims=Dict(q1=>31,coupler=>31,q2=>31),
    states=Dict(q1=>0:3,coupler=>0:4,q2=>0:3),
    parameters=Dict(ϕ1=>0.05,ϕc=>0.20,ϕ2=>0.08,EC1c=>26.0,EC2c=>26.0,EC12=>500.0))
n1, nc, n2 = recipe.inputs.q1.operators.charge, recipe.inputs.coupler.operators.charge, recipe.inputs.q2.operators.charge
J1c = 8EC1*ECc/EC1c
J2c = 8EC2*ECc/EC2c
J12 = 8EC1*EC2*(1/EC12+ECc/(EC1c*EC2c))
add_interaction!(recipe,Interaction(:q1_coupler,J1c*n1*nc))
add_interaction!(recipe,Interaction(:q2_coupler,J2c*n2*nc))
add_interaction!(recipe,Interaction(:q1_q2,J12*n1*n2))

# A full-value flux excursion: offset 0.20 is the idle operating point.
flux_pulse = Pulse(GaussianZero(4.0;duration=20);duration=20,amplitude=0.08,offset=0.20,delay=2)
flux_control = Control(:coupler_flux,ϕc,flux_pulse)
flux_gate = build(device,GateRecipe(recipe;controls=[flux_control],duration=24))
H_flux = make_hamiltonian(flux_gate)
H_flux(nothing,10.0)
check_endpoints(flux_gate;atol=1e-14) # Explicit, optional reference-endpoint check.

# Charge is exposed by the input even though a drive amplitude is absent from H.
# Extra terms are projected in the same reference basis, without rediagonalizing.
Ω = Symbolics.variable(:microwave_amplitude)
drive_recipe = draft(recipe)
drive_recipe.parameters[Ω] = 0.0
drive_term = Interaction(:q1_drive,Ω*drive_recipe.inputs.q1.operators.charge)
microwave_pulse = Pulse(GaussianZero(4.0;duration=20),4.8;duration=20,amplitude=0.03,delay=2)
drive_control = Control(:q1_microwave,Ω,microwave_pulse)
microwave_gate = build(device,GateRecipe(drive_recipe;controls=[drive_control],terms=[drive_term],duration=24))
H_microwave = make_hamiltonian(microwave_gate)
H_microwave(nothing,10.0)

both = build(device,GateRecipe(drive_recipe;controls=[flux_control,drive_control],terms=[drive_term],duration=24))
H_both = make_hamiltonian(both)
H_both(nothing,10.0)
println("Prepared dimension: ",length(both.model.spectrum),"; pulse stops: ",pulse_tstops(both))

# With no free_parameters keyword, runtime p is ignored and controls use captured
# pulses. Build parameters override recipe assignments and component idle values.
# draft(both.recipe) returns independent editable model/control/term provenance.
# Preparation captures physical definitions for calibration provenance.
# Arbitrary callable captured state remains outside the snapshot boundary.

# Evolve only the computational inputs through the full retained Hamiltonian.
# Label order follows the endpoint reference model's physical component order.
computational = [(0,0,0),(1,0,0),(0,0,1),(1,0,1)]
evolution = evolve(both;states=computational,abstol=1e-8,reltol=1e-8)
physical_matrix = gate_matrix(evolution)
reference_phases = idle_phases(evolution) # unwrapped 2π E T, label callable
println("Computational matrix size: ",size(physical_matrix),
    "; ground idle phase: ",reference_phases((0,0,0)))
# Projected columns are never normalized: their lost norm measures leakage.
# For a full unitary use evolve(both), then unitary(result).
# Histories are opt-in: evolve(both;states=computational,saveat=0:2:24).


# A small analytic X-gate calibration demonstrates optimizer choice and explicit
# acceptance without repeatedly optimizing the larger three-component model.
using OptimizationOptimJL
calibration_device = QuantumDevice(:calibration_example)
calibration_qubit = add_component!(calibration_device,make_qubit(:q;frequency=0.0))
calibration_reference = Recipe(calibration_device,[calibration_qubit])
drive_amplitude = Symbolics.variable(:drive_amplitude)
calibration_model = draft(calibration_reference)
calibration_model.parameters[drive_amplitude] = 0.0
calibration_term = Interaction(:drive,
    drive_amplitude*calibration_model.inputs.q.operators.x)
working_gate = QuantumGate(:X,GateRecipe(calibration_model;
    controls=[Control(:drive,drive_amplitude,Pulse(Constant(;duration=1.0);duration=1.0,amplitude=0.1))],
    terms=[calibration_term]))
prepared_calibration = build(calibration_device,working_gate.recipe)
x_target = ComplexF64[0 1; 1 0]
calibration = CalibrationProblem(prepared_calibration,["drive/amplitude"],x_target;
    states=[(0,),(1,)],
    evolution_settings=(;abstol=1e-10,reltol=1e-10))
run_calibrations = get(ENV,"QUANTUMDEVICES_RUN_CALIBRATIONS","false") == "true"
if run_calibrations
solution = solve(calibration,NelderMead();maxiters=100)
candidate_recipe = solution.recipe
accepted = validate_calibration(solution;save_evolution=true)
working_gate = QuantumGate(:X,accepted)
accepted_evolution = evolve(working_gate.calibration)
println("Accepted calibration loss: ",accepted.validation.loss)

end # run_calibrations

# Optimization coordinates live on calibration; the gate stays physical.
Symbolics.@variables pulse_amplitude pulse_duration
linked_calibration = CalibrationProblem(prepared_calibration,x_target;
    variables=[pulse_amplitude=>0.1,pulse_duration=>1.0],
    links=Dict("drive/amplitude"=>pulse_amplitude,"drive/duration"=>pulse_duration),
    lb=[0.0,0.5],ub=[0.5,3.0],states=[(0,),(1,)])
longer_definition = calibrated_recipe(linked_calibration,[0.125,2.0])
@assert longer_definition.controls[1].pulse.duration == 2.0
@assert isequal(linked_calibration.coordinates,(pulse_amplitude,pulse_duration))
@assert !hasproperty(longer_definition,:links)
