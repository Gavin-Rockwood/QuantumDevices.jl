# # Transmon resonator control
#
# This tutorial constructs a coupled transmon–resonator model, calibrates transmon
# π pulses, finds a sideband's drive-induced Stark shift, and calibrates its
# carrier phase. The three pulses connect the dressed states along the path
#
# $$|g,0\rangle \xrightarrow{\pi_{ge}} |e,0\rangle
#   \xrightarrow{\pi_{ef}} |f,0\rangle
#   \xrightarrow{\mathrm{sideband}} |g,1\rangle.$$
#
# Transmon rotations and charge-driven sidebands are building blocks of the
# control protocols in Huang et al., [*Fast Sideband Control of a Weakly Coupled
# Multimode Bosonic Memory*](https://arxiv.org/abs/2503.10623). Section III
# describes the sideband interaction; Supplement S7 explains how Floquet avoided
# crossings determine its resonance and rate, and how pulse ramps affect leakage.
# Here we work with one cavity mode and sin² ramps.
#
# ## 1. Set up the environment
#
# Run the code blocks in order, or run this source file from the repository root:
#
# ```sh
# julia --project=docs -e 'using Pkg; Pkg.instantiate()'
# julia --project=docs docs/src/tutorials/transmon_resonator_control.jl
# ```
#
# The calibrations and frequency sweeps can take several minutes. Time is in
# **ns**, frequencies and Hamiltonian coefficients are in **GHz**, and carrier
# phases are in radians. The model-aware `sesolve`, projected `numerical`, and
# `find_resonance` methods apply the `2π` conversion for evolution internally.

using QuantumDevices, QuantumToolbox
using OptimizationOptimJL

# ## 2. Build the coupled model
#
# The transmon Hamiltonian is $H_t=4E_C(\hat n-n_g)^2-E_J\cos\hat\varphi$.
# We couple it to a resonator with $H_r=\nu_r a^\dagger a$ through
# $H_{\mathrm{int}}=g\hat n\,i(a-a^\dagger)$.
# `op` names component operators and `param` names a model parameter.

EC = 0.10283303447280807
EJ = 26.96976142643705
resonator_frequency = 6.228083962082612
coupling = 0.026877206812551357
charge_cutoff = 60
transmon_levels = 10
resonator_levels = 10

transmon = make_transmon("transmon", EC, EJ, 2charge_cutoff + 1; ng=0.0)
resonator = make_resonator("resonator", resonator_frequency, resonator_levels)
interaction = param(:g) * op(:transmon_charge) *
    (1im * (op(:resonator_a) - op(:resonator_adag)))
model = make_model([transmon, resonator], interaction, (; g=coupling);
    truncation_dimensions=Dict(transmon => transmon_levels, resonator => resonator_levels),
    max_dimension=transmon_levels * resonator_levels)
size(model.H)

# We retain ten transmon levels and ten resonator levels, giving a 100-dimensional
# model. `model.states` contains dressed eigenstates of the coupled device,
# labeled by their connection to uncoupled transmon and resonator excitations.
# `model.spectrum` stores their energies in GHz, including static coupling shifts.
#
# | Label | Dressed state |
# |:------|:--------------|
# | `(0,0)` | $\lvert g,0\rangle$ |
# | `(1,0)` | $\lvert e,0\rangle$ |
# | `(2,0)` | $\lvert f,0\rangle$ |
# | `(0,1)` | $\lvert g,1\rangle$ |

ge_states = [model.states[(0, 0)], model.states[(1, 0)]]
ef_states = [model.states[(1, 0)], model.states[(2, 0)]]
sideband_states = [model.states[(2, 0)], model.states[(0, 1)]]

ge_frequency = abs(model.spectrum[(1, 0)] - model.spectrum[(0, 0)])
ef_frequency = abs(model.spectrum[(2, 0)] - model.spectrum[(1, 0)])

# ## 3. Calibrate the transmon π pulses
#
# ### Drive the $g\leftrightarrow e$ transition
#
# A charge drive adds $A\,s(t)\sin(2\pi\nu t+\phi)\hat n$ to the Hamiltonian.
# `Pulse` combines the amplitude, envelope, carrier, and duration. The Gaussian
# has unit peak and is centered at half the duration; its width `sigma` is in ns.
# The amplitude below is a starting guess. The gate inherits the pulse duration.

qubit_duration = 92.96875
qubit_sigma = qubit_duration / 4
ge_drive = Pulse(Gaussian(qubit_sigma);
    amplitude=0.004573, duration=qubit_duration,
    carrier=SineCarrier(ge_frequency; phase=0.0))
model.gates[:t_ge_0] = DeviceGate((; drive=ge_drive),
    param(:drive) * op(:transmon_charge))

# Start in $|g,0\rangle$, evolve the full model, and measure the amplitude in
# $|e,0\rangle$. Its squared magnitude is the target-state population.
# Keep this trial trajectory to compare with the calibrated pulse later.

evolution_options = (; abstol=1e-8, reltol=1e-8)
qubit_times = range(0, qubit_duration; length=201)
result = sesolve(model, model.gates[:t_ge_0], ge_states[1], qubit_times;
    progress_bar=false, evolution_options...)
amplitudes = state_amplitudes(Dict(:target => ge_states[2]), result)
ge_before = abs2.(amplitudes[:target])

# ### Compare with an X gate
#
# A population trace tests one input state. An X gate must also act correctly on
# superpositions, including their relative phases. Supply both ordered basis
# states to `numerical` to evolve and then project onto this two-state subspace.
# The result is a 2×2 operator; use `.data` to inspect its matrix.
#
# Use `unitary_fidelity(sigmax(), U)` with an unscaled target. For two states it
# computes $|\operatorname{Tr}(X^\dagger U)|^2/4$, ignoring only a common global
# phase. `QuantumToolbox.fidelity` measures state fidelity and should not be
# applied directly to gate matrices.
#
# We use `frame=:lab` consistently in calibration and evaluation, so the pulse
# phase also compensates relative idle evolution. Use `frame=:interaction` in
# both calls if you want to remove idle evolution from the comparison.

U_ge = numerical(model, model.gates[:t_ge_0], ge_states;
    frame=:lab, evolution_options...)
unitary_fidelity(sigmax(), U_ge)

# ### Optimize the pulse
#
# Select amplitude, phase, and frequency by their parameter paths. Duration and
# Gaussian width stay fixed. `CalibrationProblem` defines the gate loss;
# `solve` runs the optimizer, and `calibrated_gate` reconstructs the pulse.
# Nelder–Mead uses a small initial simplex appropriate to these starting guesses.

selectors = ["drive/amplitude", "drive/carrier/phase", "drive/carrier/frequency"]
simplex = OptimizationOptimJL.Optim.AffineSimplexer(a=0.0005, b=0.0005)
algorithm = NelderMead(initial_simplex=simplex)

problem = CalibrationProblem(model, model.gates[:t_ge_0], selectors, sigmax();
    states=ge_states, frame=:lab, evolution_kwargs=evolution_options)
solution = solve(problem, algorithm; maxiters=300, g_abstol=1e-8)
model.gates[:t_ge_0] = calibrated_gate(problem, solution.u)

U_ge = numerical(model, model.gates[:t_ge_0], ge_states;
    frame=:lab, evolution_options...)
unitary_fidelity(sigmax(), U_ge)

# Check the returned fidelity and `solution.retcode`, then compare the target
# populations. The calibrated coherent fidelity is approximately 0.99994.

result = sesolve(model, model.gates[:t_ge_0], ge_states[1], qubit_times;
    progress_bar=false, evolution_options...)
amplitudes = state_amplitudes(Dict(:target => ge_states[2]), result)
ge_after = abs2.(amplitudes[:target])

styles = Dict(:trial => (; color=:gray, label="Trial"),
    :calibrated => (; color=:purple, label="Calibrated"))
ge_figure, ax = plot_trajectories(
    Dict(:trial => ge_before, :calibrated => ge_after), qubit_times, styles;
    axis=(; xlabel="Time (ns)", ylabel="Target population", title="Transmon g ↔ e"),
    legend_options=(; position=:lt))
ge_figure

# ### Repeat for the $e\leftrightarrow f$ transition
#
# Change the dressed pair, frequency, and starting amplitude. The charge operator,
# Gaussian envelope, target X, and calibration procedure are the same.

ef_drive = Pulse(Gaussian(qubit_sigma);
    amplitude=0.00385, duration=qubit_duration,
    carrier=SineCarrier(ef_frequency; phase=0.0))
model.gates[:t_ef_0] = DeviceGate((; drive=ef_drive),
    param(:drive) * op(:transmon_charge))

result = sesolve(model, model.gates[:t_ef_0], ef_states[1], qubit_times;
    progress_bar=false, evolution_options...)
amplitudes = state_amplitudes(Dict(:target => ef_states[2]), result)
ef_before = abs2.(amplitudes[:target])

problem = CalibrationProblem(model, model.gates[:t_ef_0], selectors, sigmax();
    states=ef_states, frame=:lab, evolution_kwargs=evolution_options)
solution = solve(problem, algorithm; maxiters=300, g_abstol=1e-8)
model.gates[:t_ef_0] = calibrated_gate(problem, solution.u)

U_ef = numerical(model, model.gates[:t_ef_0], ef_states;
    frame=:lab, evolution_options...)
unitary_fidelity(sigmax(), U_ef)

result = sesolve(model, model.gates[:t_ef_0], ef_states[1], qubit_times;
    progress_bar=false, evolution_options...)
amplitudes = state_amplitudes(Dict(:target => ef_states[2]), result)
ef_after = abs2.(amplitudes[:target])

ef_figure, ax = plot_trajectories(
    Dict(:trial => ef_before, :calibrated => ef_after), qubit_times, styles;
    axis=(; xlabel="Time (ns)", ylabel="Target population", title="Transmon e ↔ f"),
    legend_options=(; position=:lt))
ef_figure

# Population exchange alone does not establish an X gate: the relative phase
# must also be correct. If amplitude and frequency already give the desired
# exchange, calibration can be restricted to `["drive/carrier/phase"]`.
# We use that simpler approach for the sideband below.
#
# ## 4. Find the sideband resonance and Stark shift
#
# The charge-driven $|f,0\rangle\leftrightarrow|g,1\rangle$ sideband transfers a
# transmon excitation into a cavity photon. At strong drive, its resonance shifts
# from the undriven dressed energy difference.
#
# `find_resonance` tracks the two relevant Floquet modes and fits the minimum of
# their quasienergy separation. The minimum gives the driven resonance; the gap
# gives the constant-drive swap time, $t_\pi=1/(2\Delta_{\min})$, when the
# interaction is well described by an isolated pair. This method is described in
# [Supplement S7](https://arxiv.org/html/2503.10623v1) of the paper.
#
# First bracket the avoided crossing. Supply the full idle Hamiltonian, the
# charge operator, and the plateau amplitude. This overload constructs a
# continuous sine drive; the finite pulse envelope is omitted because Floquet
# analysis requires a periodic Hamiltonian.

sideband_amplitude = 0.735
sideband_bare_frequency = abs(model.spectrum[(2, 0)] - model.spectrum[(0, 1)])
sideband_charge = numerical(model, op(:transmon_charge))
references = Dict(:f0 => sideband_states[1], :g1 => sideband_states[2])

frequencies = sideband_bare_frequency .+ range(-0.08, 0.01; length=31)
resonance = find_resonance(model.H, sideband_charge, sideband_amplitude,
    frequencies, references; state_keys=[:f0, :g1],
    propagator_kwargs=(; abstol=1e-10, reltol=1e-10))

# This window spans −80 to +10 MHz relative to the undriven transition. It must
# contain the avoided-crossing minimum; adjust it for another device or amplitude.
# Refine the estimate in a narrower ±5 MHz window:

frequencies = range(resonance.frequency - 0.005, resonance.frequency + 0.005; length=31)
resonance = find_resonance(model.H, sideband_charge, sideband_amplitude,
    frequencies, references; state_keys=[:f0, :g1],
    propagator_kwargs=(; abstol=1e-10, reltol=1e-10))
resonance.frequency

# The signed Stark shift is the driven resonance minus the undriven transition
# frequency. Multiply GHz by 1000 to express it in MHz. Here it is approximately
# **−37.654 MHz**, with a driven resonance near **2.832308 GHz**.

stark_shift = resonance.frequency - sideband_bare_frequency
1000 * stark_shift

# Plot the spectrum and avoided crossing to check that the fit follows the data:

resonance_figure, resonance_axes = plot_resonance(resonance;
    frequency_offset=sideband_bare_frequency,
    xlabel="Drive detuning from undriven dressed transition (GHz)")
resonance_figure

# ## 5. Calibrate the sideband SWAP phase
#
# Use the fitted frequency with a 205.054 ns pulse and 10 ns sin² ramps. The
# constant-drive swap time is a starting estimate: ramps change the effective
# interaction time, so the complete pulse still needs a time-domain check.
# During the drive, population can temporarily leave the two undriven reference
# states as they become dressed. The ramp down can return that population.

sideband_duration = 205.05411739244443
sb_drive = Pulse(RampedFlattop(10.0; ramp=:sine_squared);
    amplitude=sideband_amplitude, duration=sideband_duration,
    carrier=SineCarrier(resonance.frequency; phase=0.0))
model.gates[:sb_f0g1] = DeviceGate((; drive=sb_drive),
    param(:drive) * op(:transmon_charge))

sideband_times = range(0, sideband_duration; length=1001)
result = sesolve(model, model.gates[:sb_f0g1], sideband_states[1], sideband_times;
    progress_bar=false, evolution_options...)
amplitudes = state_amplitudes(Dict(:target => sideband_states[2]), result)
sb_before = abs2.(amplitudes[:target])

U_sb = numerical(model, model.gates[:sb_f0g1], sideband_states;
    frame=:lab, evolution_options...)
unitary_fidelity(sigmax(), U_sb)

# Population exchange is already good, but the coherent X overlap is only about
# 0.4135. Optimize only the carrier phase, holding the other pulse parameters fixed:

problem = CalibrationProblem(model, model.gates[:sb_f0g1],
    ["drive/carrier/phase"], sigmax();
    states=sideband_states, frame=:lab, evolution_kwargs=evolution_options)
solution = solve(problem, algorithm; maxiters=300, g_abstol=1e-8)
model.gates[:sb_f0g1] = calibrated_gate(problem, solution.u)

U_sb = numerical(model, model.gates[:sb_f0g1], sideband_states;
    frame=:lab, evolution_options...)
unitary_fidelity(sigmax(), U_sb)

# The coherent fidelity improves to approximately 0.9996. The two population
# curves remain very similar: this improvement primarily corrects their phases.

result = sesolve(model, model.gates[:sb_f0g1], sideband_states[1], sideband_times;
    progress_bar=false, evolution_options...)
amplitudes = state_amplitudes(Dict(:target => sideband_states[2]), result)
sb_after = abs2.(amplitudes[:target])

sb_figure, ax = plot_trajectories(
    Dict(:trial => sb_before, :calibrated => sb_after), sideband_times, styles;
    axis=(; xlabel="Time (ns)", ylabel="Target population", title="Sideband f,0 ↔ g,1"),
    legend_options=(; position=:lt))
sb_figure

# The target X acts on `[|f,0⟩, |g,1⟩]`, rather than the entire device Hilbert
# space. Projection preserves leakage, so the resulting 2×2 matrix need not be
# exactly unitary. These closed-system gate overlaps exclude decoherence and
# experimental preparation and readout errors.
#
# The calibrated pulses are stored in `model.gates`. See [Persistence](../user_guide/persistence.md)
# to save and reload the model for further work.
