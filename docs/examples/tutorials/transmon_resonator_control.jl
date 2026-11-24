# # Transmon–resonator control
#
# We build a coupled transmon and resonator, calibrate transmon population transfer,
# and use a driven sideband to move an excitation into the resonator:
#
# $$|g,0\rangle\xrightarrow{\pi_{ge}}|e,0\rangle
# \xrightarrow{\pi_{ef}}|f,0\rangle\xrightarrow{\mathrm{sideband}}|g,1\rangle.$$
#
# This is the original tutorial's physical setup, migrated to Recipe, Control
# and the current physical gate evolution. Charge-driven sidebands are discussed
# in Huang et al., [Fast Sideband Control of a Weakly Coupled Multimode Bosonic
# Memory](https://arxiv.org/abs/2503.10623), Section III and Supplement S7.
# Normal docs builds execute only the waveform plot. To include calibration,
# population and Floquet figures, build with `DOCS_RUN_TUTORIALS=true` (see below).
# Historical fidelities are not reproduced as new results.
#
# ## 1. Environment and execution
#
# Time is ns, energies and frequencies are GHz (`E/h`), and phase is radians.
# Built-in `evolve` and Floquet helpers apply the `2π` conversion once.
# To run the full numerical workflow explicitly:
#
# ```sh
# julia --project=docs docs/examples/tutorials/transmon_resonator_control.jl
# ```
#
# Invoking the source deliberately runs the full 100-dimensional workflow,
# including long trajectories, scans and calibrations. Documentation builds render
# this code without executing it; only the small waveform example runs by default.
# To execute the complete tutorials and capture their figures in the website:
#
# ```sh
# DOCS_RUN_TUTORIALS=true julia --startup-file=no --project=docs docs/make.jl
# ```
#
# Stop `docs:dev` first: draft mode skips execution and cannot use this flag.
# After the build, use `npm --prefix docs run docs:preview` to view the results.

using QuantumDevices, QuantumToolbox, CairoMakie
using OptimizationOptimJL
import Symbolics

# ## 2. Define the coupled device
#
# The transmon has $H_t=4E_C(n-n_g)^2-E_J\cos\varphi$, the resonator has
# $H_r=\nu_r a^\dagger a$, and $H_{int}=g n\,i(a-a^\dagger)$.
# The quadrature and sign of this interaction are part of the physical setup.

EC = 0.10283303447280807
EJ = 26.96976142643705
resonator_frequency = 6.228083962082612
coupling = 0.026877206812551357
charge_cutoff = 60
transmon_levels = 10
resonator_levels = 10
device = QuantumDevice(:transmon_resonator)
transmon = add_component!(device,make_transmon(:transmon,EC,EJ;dimension=2charge_cutoff+1,ng=0.0))
resonator = add_component!(device,make_resonator(:resonator;frequency=resonator_frequency,dimension=resonator_levels))
recipe = Recipe(device,[transmon,resonator];
    states=Dict(transmon=>0:transmon_levels-1,resonator=>0:resonator_levels-1))
n = recipe.inputs.transmon.operators.charge
a = recipe.inputs.resonator.operators.a
g = Symbolics.variable(:charge_coupling)
recipe.parameters[g] = coupling
add_interaction!(recipe,Interaction(:coupling,g*n*(1im*(a-a'))))
model = build(device,recipe)
@assert length(model.spectrum) == 100
ge_states = [(0,0),(1,0)]
ef_states = [(1,0),(2,0)]
sideband_states = [(2,0),(0,1)]
ge_frequency = abs(stateenergy(model,ge_states[2])-stateenergy(model,ge_states[1]))
ef_frequency = abs(stateenergy(model,ef_states[2])-stateenergy(model,ef_states[1]))

# Labels follow `[transmon,resonator]`: `(0,0)` is g,0, `(1,0)` is e,0,
# `(2,0)` is f,0 and `(0,1)` is g,1. Inspect label diagnostics and converge
# both the charge cutoff and retained levels before trusting quantitative results.
# A new `build(...;parameters=...)` recomputes the reference basis; a control
# varies coefficients in that fixed basis during one evolution.
#
# ## 3. Prepare the transmon π pulses
#
# A charge drive adds a term; its gate-owned recipe supplies the zero idle value.
# Gaussian remains the original truncated Gaussian, not a zero-ended replacement.
# A positional carrier uses cosine: `phase=-pi/2` preserves the original sine.

qubit_duration = 92.96875
qubit_sigma = qubit_duration/4
amplitude = Symbolics.variable(:microwave_amplitude)
drive_recipe = draft(recipe)
drive_recipe.parameters[amplitude] = 0.0
drive = Interaction(:charge_drive,amplitude*drive_recipe.inputs.transmon.operators.charge)
shape = Gaussian(;duration=qubit_duration,sigma=qubit_sigma)
pulse_times = range(0,qubit_duration;length=201)
control_figure = Figure(size=(720,300))
control_axis = Axis(control_figure[1,1];xlabel="Time (ns)",ylabel="Drive envelope (GHz)")
lines!(control_axis,pulse_times,0.004573 .* shape.(pulse_times,qubit_duration);label="g ↔ e trial")
lines!(control_axis,pulse_times,0.00385 .* shape.(pulse_times,qubit_duration);label="e ↔ f trial")
axislegend(control_axis; position=:lt)
control_figure

# cheap-docs-begin
# ```@example transmon_controls
# using QuantumDevices, CairoMakie
# T = 92.96875
# e = Gaussian(;duration=T,sigma=T/4)
# t = range(0,T;length=201)
# fig = Figure(size=(720,300))
# ax = Axis(fig[1,1];xlabel="Time (ns)",ylabel="Drive envelope (GHz)")
# lines!(ax,t,.004573 .* e.(t,T);label="g ↔ e trial")
# lines!(ax,t,.00385 .* e.(t,T);label="e ↔ f trial")
# axislegend(ax; position=:lt)
# fig
# ```
# cheap-docs-end
#
# ### Evolve both inputs and plot a population
#
# `states` chooses ordered input columns, not a reduced dynamical Hamiltonian.
# Both columns evolve in the full retained space. `gate_matrix` projects them
# without renormalizing, preserving leakage and physical relative phases.

population_exchange_cost(target,M) = 1-(abs2(M.data[2,1])+abs2(M.data[1,2]))/2
x_target = ComplexF64[0 1;1 0]
evolution_options = (;abstol=1e-8,reltol=1e-8,energy_shift=:auto)
ge_pulse = Pulse(shape,ge_frequency;amplitude=0.004573,phase=-pi/2)
ge_recipe = GateRecipe(drive_recipe;terms=[drive],controls=[Control(:drive,amplitude,ge_pulse)])
ge_prepared = build(device,ge_recipe)
ge_trial = evolve(ge_prepared;states=ge_states,saveat=pulse_times,evolution_options...)
ge_before = [abs2(gate_matrix(ge_trial;step=i).data[2,1]) for i in eachindex(ge_trial.times)]
ge_trial_loss = population_exchange_cost(x_target,gate_matrix(ge_trial))

# ### Calibrate the amplitude for population exchange
#
# The objective averages transfer over both inputs without renormalizing away
# leakage. It ignores relative phases, so high transfer does not establish a
# coherent X gate. Keep phase, frequency, duration and Gaussian width fixed.
# No interaction-frame inverse is applied.

selectors = ["drive/amplitude"]
simplex = OptimizationOptimJL.Optim.AffineSimplexer(a=0.0005,b=0.0005)
algorithm = NelderMead(initial_simplex=simplex)
ge_problem = CalibrationProblem(ge_prepared,selectors,x_target;
    states=ge_states,objective=population_exchange_cost,evolution_settings=evolution_options)
ge_fit = solve(ge_problem,algorithm;maxiters=300,g_abstol=1e-8)
ge_checked = validate_calibration(ge_fit;save_evolution=true)
ge_checked.solution.retcode
ge_checked.validation.loss
ge_after_history = evolve(ge_checked;saveat=pulse_times)
ge_after = [abs2(gate_matrix(ge_after_history;step=i).data[2,1]) for i in eachindex(ge_after_history.times)]
ge_figure = Figure(size=(720,300))
ax = Axis(ge_figure[1,1];xlabel="Time (ns)",ylabel="Target population",title="Transmon g ↔ e")
lines!(ax,ge_trial.times,ge_before;label="Trial",color=:gray)
lines!(ax,ge_after_history.times,ge_after;label="Calibrated",color=:purple)
axislegend(ax; position=:lt)
ge_figure

# Inspect return codes, validation loss and leakage before adopting
# `QuantumGate(:ge,ge_checked)`. A successful optimizer return alone does not
# establish physical convergence or cutoff convergence.
#
# ### Repeat for e ↔ f
#
# Only the dressed pair, frequency and trial amplitude change. The charge term,
# envelope, target matrix and selected paths are reused.

ef_pulse = Pulse(shape,ef_frequency;amplitude=0.00385,phase=-pi/2)
ef_prepared = build(device,GateRecipe(drive_recipe;terms=[drive],controls=[Control(:drive,amplitude,ef_pulse)]))
ef_trial = evolve(ef_prepared;states=ef_states,saveat=pulse_times,evolution_options...)
ef_before = [abs2(gate_matrix(ef_trial;step=i).data[2,1]) for i in eachindex(ef_trial.times)]
ef_fit = solve(CalibrationProblem(ef_prepared,selectors,x_target;
    states=ef_states,objective=population_exchange_cost,evolution_settings=evolution_options),algorithm;maxiters=300,g_abstol=1e-8)
ef_checked = validate_calibration(ef_fit)
ef_after_history = evolve(ef_checked;saveat=pulse_times)
ef_after = [abs2(gate_matrix(ef_after_history;step=i).data[2,1]) for i in eachindex(ef_after_history.times)]
ef_figure = Figure(size=(720,300))
ax = Axis(ef_figure[1,1];xlabel="Time (ns)",ylabel="Target population",title="Transmon e ↔ f")
lines!(ax,ef_trial.times,ef_before;label="Trial",color=:gray)
lines!(ax,ef_after_history.times,ef_after;label="Calibrated",color=:purple)
axislegend(ax; position=:lt)
ef_figure

# ## 4. Locate the sideband and its Stark shift
#
# A periodic constant-amplitude drive isolates the Floquet avoided crossing;
# the finite Gaussian pulse is not a periodic drive. Start with the original
# bracket, −80 to +10 MHz around the undriven dressed transition, then refine
# within ±5 MHz of the fitted minimum. Inspect residuals, mode tracking and
# competing crossings before interpreting the fitted frequency or gap.

sideband_amplitude = 0.735
bare_frequency = abs(stateenergy(model,sideband_states[2])-stateenergy(model,sideband_states[1]))
frequencies = bare_frequency .+ range(-0.08,0.01;length=31)
resonance = find_resonance(model,transmon.operators["charge"],sideband_amplitude,frequencies;
    states=sideband_states,phase=-pi/2,evolution_settings=(;abstol=1e-10,reltol=1e-10))
refined_frequencies = range(resonance.frequency-0.005,resonance.frequency+0.005;length=31)
resonance = find_resonance(model,transmon.operators["charge"],sideband_amplitude,refined_frequencies;
    states=sideband_states,phase=-pi/2,evolution_settings=(;abstol=1e-10,reltol=1e-10))
shift = stark_shift(resonance)
shift
resonance_figure,resonance_axes = plot_resonance(resonance;frequency_offset=bare_frequency)
resonance_figure

# The signed shift is driven resonance minus undriven transition frequency.
# Multiply GHz by 1000 for MHz. The approximate constant-drive swap time
# `1/(2*minimum_gap)` is a starting estimate; finite ramps change the effective
# interaction time. No historical numerical shift is asserted here.
#
# ## 5. Calibrate the finite sideband plateau time
#
# Start from the original 205.054 ns duration and retain 10 ns sin² ramps.
# Fit only the plateau time for population exchange, keeping amplitude, carrier
# frequency and phase fixed. Population may leave the undriven pair during the
# plateau and return during ramp-down. This fit does not constrain relative phases.

sideband_duration = 205.05411739244443
sideband_shape = RampedFlattop(10.;duration=sideband_duration,ramp=:sine_squared)
sb_pulse = Pulse(sideband_shape,resonance.frequency;amplitude=sideband_amplitude,phase=-pi/2)
sb_prepared = build(device,GateRecipe(drive_recipe;terms=[drive],controls=[Control(:drive,amplitude,sb_pulse)]))
sb_times = range(0,sideband_duration;length=1001)
sb_trial = evolve(sb_prepared;states=sideband_states,saveat=sb_times,evolution_options...)
sb_before = [abs2(gate_matrix(sb_trial;step=i).data[2,1]) for i in eachindex(sb_trial.times)]
sb_plateau_time = Symbolics.variable(:sb_plateau_time)
sb_links = Dict("drive/duration"=>sb_plateau_time+20.,
    "drive/envelope/rise_time"=>10.,"drive/envelope/fall_time"=>10.)
sb_fit = solve(CalibrationProblem(sb_prepared,x_target;
    variables=[sb_plateau_time=>sideband_duration-20.],links=sb_links,
    states=sideband_states,objective=population_exchange_cost,evolution_settings=evolution_options),algorithm;maxiters=300,g_abstol=1e-8)
sb_checked = validate_calibration(sb_fit)
sb_after_times = range(0,sb_checked.recipe.duration;length=1001)
sb_after_history = evolve(sb_checked;saveat=sb_after_times)
sb_after = [abs2(gate_matrix(sb_after_history;step=i).data[2,1]) for i in eachindex(sb_after_history.times)]
sb_figure = Figure(size=(720,300))
ax = Axis(sb_figure[1,1];xlabel="Time (ns)",ylabel="Target population",title="Sideband f,0 ↔ g,1")
lines!(ax,sb_trial.times,sb_before;label="Trial",color=:gray)
lines!(ax,sb_after_history.times,sb_after;label="Calibrated",color=:purple)
axislegend(ax; position=:lt)
sb_figure

# Capture the accepted sideband gate after inspecting the population plot.

sideband_gate = QuantumGate(:sideband,sb_checked)

# These targets act on ordered dressed pairs, not the entire Hilbert space.
# Report projected loss and missing column norm separately. The closed-system
# model omits decoherence and experimental preparation/readout errors.
# Keep accepted QuantumGate results with the device using the
# [persistence workflow](../user_guide/persistence.md); captured calibration
# settings and physical context make later replay explicit.
