# # Tunable coupler control
#
# Two flux-tunable transmons exchange an excitation through a third transmon:
#
# $$|1,0,0\rangle\longleftrightarrow|0,1,0\rangle.$$
#
# We plot the flux controls, simulate transfer, fit population exchange and then
# inspect the full computational subspace. A high exchange probability does not
# establish a coherent SWAP: relative phases and the two-excitation sector matter.
# This preserves the original tutorial's physical values and tensor order.
# Old reported scores are not presented as measurements from this API.
#
# ## 1. Environment
#
# Run this full source deliberately: it builds the original model and performs
# calibrations, which can take several minutes. A documentation build renders
# the code but executes only the small waveform plot below by default.
# To capture all population and calibration figures in the website, stop
# `docs:dev` and run:
#
# ```sh
# DOCS_RUN_TUTORIALS=true julia --startup-file=no --project=docs docs/make.jl
# npm --prefix docs run docs:preview
# ```
#
# Draft mode skips execution and cannot be combined with this flag.
#
# ```sh
# julia --project=docs docs/examples/tutorials/tunable_coupler_control.jl
# ```
#
# Time is ns, energies are GHz (`E/h`), flux is Φext/Φ₀, and phases are radians.
# `evolve` applies `2π` once; do not multiply its Hamiltonian a second time.

using QuantumDevices, QuantumToolbox, CairoMakie
using OptimizationOptimJL
import Symbolics

# ## 2. Build the coupled model
#
# Each symmetric split transmon has effective Josephson energy
# $E_J(\phi)=E_{J,max}|\cos(\pi\phi)|$. The current constructor accepts the
# total EJmax and asymmetry `d=0`, rather than two equal junction energies.

EC1, EC2, ECc = 0.198, 0.18, 0.097
EJ1, EJ2, EJc = 11.143, 11.15, 33.528
EC1c, EC2c, EC12 = 4.3, 4.5, 161.4
charge_cutoff = 60
retained_levels = 3
device = QuantumDevice(:tunable_coupler)
q1 = add_component!(device,make_flux_tunable_transmon(:q1,EC1,EJ1;dimension=2charge_cutoff+1,d=0.,flux=0.))
q2 = add_component!(device,make_flux_tunable_transmon(:q2,EC2,EJ2;dimension=2charge_cutoff+1,d=0.,flux=0.))
coupler = add_component!(device,make_flux_tunable_transmon(:c,ECc,EJc;dimension=2charge_cutoff+1,d=0.,flux=0.))
recipe = Recipe(device,[q1,q2,coupler];states=Dict(q1=>0:retained_levels-1,
    q2=>0:retained_levels-1,coupler=>0:retained_levels-1))

# Capacitive coupling uses the same circuit parameters and signs as before.
# Interactions refer to this recipe's scoped operators. Registering parameters
# separately makes their physical identities explicit. EC1c, EC2c and EC12
# are circuit charging-energy scales, not the direct interaction coefficients:
#
# $$g_{1c}=8E_{C1}E_{Cc}/E_{C1c},\quad g_{2c}=8E_{C2}E_{Cc}/E_{C2c},$$
# $$g_{12}=8E_{C1}E_{C2}(1/E_{C12}+E_{Cc}/(E_{C1c}E_{C2c})).$$

Symbolics.@variables g1c g2c g12
η = EC12*ECc/(EC1c*EC2c)
add_parameters!(device,Dict(g1c=>8EC1*ECc/EC1c,g2c=>8EC2*ECc/EC2c,
    g12=>8(1+η)*EC1*EC2/EC12))
n1 = recipe.inputs.q1.operators.charge
n2 = recipe.inputs.q2.operators.charge
nc = recipe.inputs.c.operators.charge
add_interaction!(recipe,Interaction(:q1_coupler,g1c*n1*nc))
add_interaction!(recipe,Interaction(:q2_coupler,g2c*n2*nc))
add_interaction!(recipe,Interaction(:q1_q2,g12*n1*n2))
model = build(device,recipe)
@assert length(model.spectrum) == retained_levels^3
exchange_states = [(1,0,0),(0,1,0)]

# Each component starts in a 121-state charge basis and retains three levels.
# Labels follow `[q1,q2,coupler]`; `(1,0,0)` and `(0,1,0)` refer to dressed
# idle states connected to a single excitation in each qubit. Increase both
# cutoffs and retained levels before using the result quantitatively.
#
# ## 3. Define and inspect the flux pulse
#
# The q1 excursion lowers its frequency; the coupler excursion changes mediated
# coupling. Keep the original duration/ramp times fixed during amplitude fitting.
# The coupler starts one ramp later and finishes one ramp earlier than q1.

pulse_duration = 28.135660302912168
ramp_time = 1.8727503955158111
q1_flux = Pulse(RampedFlattop(ramp_time;duration=pulse_duration);amplitude=0.130,offset=0.)
coupler_flux = Pulse(RampedFlattop(ramp_time;duration=pulse_duration-2ramp_time);
    amplitude=0.385,offset=0.,delay=ramp_time)
flux1,fluxc = q1.variables["flux"],coupler.variables["flux"]
exchange_recipe = GateRecipe(recipe;controls=[Control(:q1_flux,flux1,q1_flux),
    Control(:coupler_flux,fluxc,coupler_flux)])
prepared = build(device,exchange_recipe)
@assert check_endpoints(prepared;atol=1e-14)
times = range(0,pulse_duration;length=201)
flux_figure = Figure(size=(720,300))
ax = Axis(flux_figure[1,1];xlabel="Time (ns)",ylabel="Flux (Φ₀)",title="Trial exchange controls")
lines!(ax,times,q1_flux.(times);label="Qubit 1",color=:royalblue)
lines!(ax,times,coupler_flux.(times);label="Coupler",color=:purple)
axislegend(ax; position=:lt)
flux_figure

# cheap-docs-begin
# The plot is also generated cheaply during documentation builds, without building
# the model or evolving the gate:
#
# ```@example coupler_controls
# using QuantumDevices, CairoMakie
# T = 28.135660302912168
# r = 1.8727503955158111
# q = Pulse(RampedFlattop(r;duration=T);amplitude=.130)
# c = Pulse(RampedFlattop(r;duration=T-2r);amplitude=.385,delay=r)
# t = range(0,T;length=201)
# fig = Figure(size=(720,300))
# ax = Axis(fig[1,1];xlabel="Time (ns)",ylabel="Flux (Φ₀)",title="Original trial controls")
# lines!(ax,t,q.(t);label="Qubit 1",color=:royalblue)
# lines!(ax,t,c.(t);label="Coupler",color=:purple)
# axislegend(ax; position=:lt)
# fig
# ```
# cheap-docs-end
#
# ## 4. Simulate transfer and preserve leakage
#
# Selected columns still evolve in the full retained model. Project each saved
# propagator onto the ordered exchange pair to obtain target population.

evolution_options = (;abstol=1e-8,reltol=1e-8)
trial = evolve(prepared;states=exchange_states,saveat=times,evolution_options...)
transfer_before = [abs2(gate_matrix(trial;step=i).data[2,1]) for i in eachindex(trial.times)]
trial_matrix = gate_matrix(trial)
trial_matrix.data

# Mean population transfer over both directions is a separate objective. Missing
# column norm is not renormalized away. This objective intentionally ignores
# relative phases; the coherent gate loss below does not.

population_loss(target,M) = 1-(abs2(M.data[2,1])+abs2(M.data[1,2]))/2
x_target = ComplexF64[0 1;1 0]
trial_transfer = 1-population_loss(x_target,trial_matrix)

# ## 5. Calibrate the two flux amplitudes
#
# Select only those amplitudes. Duration, timing and physical idle model stay
# fixed. Retain the original nearby Nelder–Mead simplex and iteration budget.

selectors = ["q1_flux/amplitude","coupler_flux/amplitude"]
simplex = OptimizationOptimJL.Optim.AffineSimplexer(a=0.001,b=0.001)
algorithm = NelderMead(initial_simplex=simplex)
problem = CalibrationProblem(prepared,selectors,x_target;states=exchange_states,
    objective=population_loss,evolution_settings=evolution_options)
fit = solve(problem,algorithm;maxiters=100,g_abstol=1e-8)
checked = validate_calibration(fit)
checked.solution.retcode
checked.validation.loss
calibrated_history = evolve(checked;saveat=times)
transfer_after = [abs2(gate_matrix(calibrated_history;step=i).data[2,1]) for i in eachindex(calibrated_history.times)]
transfer_figure = Figure(size=(720,300))
ax = Axis(transfer_figure[1,1];xlabel="Time (ns)",ylabel="Target population",title="Coupler exchange")
lines!(ax,trial.times,transfer_before;label="Trial",color=:gray)
lines!(ax,calibrated_history.times,transfer_after;label="Calibrated",color=:purple)
axislegend(ax; position=:lt)
transfer_figure

# Validate cutoff/solver convergence and inspect the return code before accepting
# the result. Optimizer convergence and physical gate quality are separate checks.
#
# ### Check the full computational action
#
# The computational order is |00⟩, |01⟩, |10⟩, |11⟩ with coupler in its ground
# state. Replay the accepted candidate with all four inputs, then compare to SWAP.

computational = [(0,0,0),(0,1,0),(1,0,0),(1,1,0)]
computational_history = evolve(checked;states=computational)
M = gate_matrix(computational_history)
swap_target = ComplexF64[1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1]
coherent_swap_loss = gate_infidelity(swap_target,M)
leakage = 1 .- vec(sum(abs2,M.data;dims=1))
(coherent_swap_loss,leakage)

# Population fitting need not produce a coherent SWAP. Relative phases and the
# action on |11⟩ remain unconstrained by the two-input population objective.
# Current evolution keeps physical lab phases; no interaction-frame inverse or
# virtual-Z correction is silently applied.
#
# ## 6. Add microwave control
#
# A charge drive adds a term rather than replacing a flux variable. Preserve the
# original π/2 starting pulse; the phase is shifted by −π/2 to express its sine
# carrier with the current cosine convention. Calibrate on both spectator states
# when targeting a single-qubit gate in the two-qubit computational space.

q1_frequency = abs(stateenergy(model,(1,0,0))-stateenergy(model,(0,0,0)))
drive_amplitude = Symbolics.variable(:microwave_amplitude)
microwave_term = Interaction(:q1_drive,drive_amplitude*n1;parameters=Dict(drive_amplitude=>0.))
microwave = Pulse(SineSquared(;duration=9.027992059968785),q1_frequency;
    amplitude=0.024775449345990953,phase=1.5501484428191823-pi/2)
microwave_recipe = GateRecipe(recipe;terms=[microwave_term],
    controls=[Control(:drive,drive_amplitude,microwave)])
microwave_prepared = build(device,microwave_recipe)

# Store an accepted result as `QuantumGate(:exchange,checked)` only after inspecting
# its validation. The [persistence guide](../user_guide/persistence.md) explains
# device-owned gates and portable captured physical context. These simulations
# omit decoherence, preparation errors and measurement errors.
