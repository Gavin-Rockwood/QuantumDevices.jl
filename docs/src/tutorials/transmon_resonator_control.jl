# # Transmon resonator control
#
# In this walkthrough we build a transmon coupled to a resonator, define three
# charge drives, and calibrate their population transfers. Run each code block in
# order, just like cells in a notebook.
#
# ## 1. Define the device and build the model
#
# From the repository root, start Julia with the docs environment:
#
# ```sh
# julia --project=docs -e 'using Pkg; Pkg.instantiate()'
# julia --project=docs
# ```
#
# Time is in ns and the frequencies below are in GHz (cycles/ns). Hamiltonian
# coefficients stay in GHz. Multiply the Hamiltonian by `2π` at evolution.

using QuantumDevices, QuantumToolbox, CairoMakie, LinearAlgebra, SciMLBase
using OptimizationOptimJL

EC = 0.10283303447280807
EJ = 26.96976142643705
resonator_frequency = 6.228083962082612
coupling = 0.026877206812551357
charge_cutoff = 60
transmon_levels = 10
resonator_levels = 10

q = make_transmon("q", EC, EJ, 2charge_cutoff + 1; ng=0.0)
r = make_resonator("r", resonator_frequency, resonator_levels)
interaction = param(:g) * op(:q_charge) * (1im * (op(:r_a) - op(:r_adag)))
model = make_model([q, r], interaction, (; g=coupling);
    truncation_dimensions=Dict(q => transmon_levels, r => resonator_levels),
    max_dimension=transmon_levels * resonator_levels)
size(model.H)

# The transmon uses a 121-state charge basis before retaining ten levels. Together
# with ten resonator levels, this gives a 100-dimensional coupled model.
# `model.states[(i, j)]` gives a dressed state: `(0, 0)` is `|g,0⟩`, `(1, 0)` is
# `|e,0⟩`, `(2, 0)` is `|f,0⟩`, and `(0, 1)` is `|g,1⟩`.
#
# ## 2. Define the three drives
#
# Use Gaussian envelopes for the qubit transitions and a flattop with sin² ramps
# for `|f,0⟩→|g,1⟩`. Each sideband ramp lasts 11.6257 ns.

qubit_duration = 92.96875
qubit_sigma = qubit_duration / 4
sideband_duration = 205.05411739244443
ramp_time = 11.6257

ge_drive = Pulse(Gaussian(qubit_sigma);
    amplitude=0.004573, duration=qubit_duration,
    carrier=SineCarrier(4.6051057859752005))
ef_drive = Pulse(Gaussian(qubit_sigma);
    amplitude=0.0032725, duration=qubit_duration,
    carrier=SineCarrier(4.496675490739561))
sb_drive = Pulse(RampedFlattop(ramp_time);
    amplitude=0.7203, duration=sideband_duration,
    carrier=SineCarrier(-2.832053845801255))

model.gates[:q_ge_0] = DeviceGate((; drive=ge_drive),
    param(:drive) * op(:q_charge))
model.gates[:q_ef_0] = DeviceGate((; drive=ef_drive),
    param(:drive) * op(:q_charge))
model.gates[:sb_f0g1] = DeviceGate((; drive=sb_drive),
    param(:drive) * op(:q_charge))

# These amplitudes and frequencies are starting guesses. Calibration will adjust
# both while keeping durations, Gaussian width, and ramp time fixed.
# Each `Pulse` owns its duration; the gate infers its duration from its controls.
# The envelope and carrier are reusable objects with discoverable parameters.
# Each drive is `amplitude × envelope × sin(2π × frequency × t)`.
#
# Inspect the sideband envelope before running the simulation:

sideband_envelope = Pulse(sb_drive.envelope; duration=sideband_duration)
envelope_times = range(0, sideband_duration; length=501)
controls = Figure(size=(700, 280))
ax_envelope = Axis(controls[1, 1]; xlabel="Time (ns)", ylabel="Envelope",
    title="Sideband envelope with sin² ramps")
lines!(ax_envelope, envelope_times,
    [sideband_envelope(t) for t in envelope_times])
controls

# ## 3. Simulate the trial pulses
#
# For each transition, evolve the initial dressed state with `sesolve` and measure
# its overlap with the target. The small function below also measures population
# outside the initial–target pair. Dividing by the evolved state's squared norm
# removes small numerical norm drift.

function populations(model, gate, initial, target; samples=201)
    times = range(0, gate.duration; length=samples)
    result = sesolve(2pi * numerical(model, gate), initial, times;
        progress_bar=false, tstops=pulse_tstops(gate), abstol=1e-8, reltol=1e-8)
    @assert SciMLBase.successful_retcode(result.retcode)
    norms = [real(dot(state, state)) for state in result.states]
    @assert maximum(abs.(norms .- 1)) < 1e-4
    initial_population = [abs2(dot(initial, state)) / n
        for (state, n) in zip(result.states, norms)]
    target_population = [abs2(dot(target, state)) / n
        for (state, n) in zip(result.states, norms)]
    leakage = max.(0.0, 1 .- initial_population .- target_population)
    return (; times, target_population, leakage)
end

transfers = [
    (:q_ge_0, (0, 0), (1, 0)),
    (:q_ef_0, (1, 0), (2, 0)),
    (:sb_f0g1, (2, 0), (0, 1)),
]
trial_gates = copy(model.gates)
before = Dict(name => populations(model, trial_gates[name],
    model.states[initial], model.states[target])
    for (name, initial, target) in transfers)
[(name, last(before[name].target_population)) for (name, _, _) in transfers]

# The trial qubit pulses transfer about 92% of the population. The trial sideband
# transfers about 59%. Next we improve those starting guesses.
#
# ## 4. Calibrate amplitude and frequency
#
# Our loss is `1 - target_population` at the end of the pulse. Select the pulse's
# amplitude and frequency with `CalibrationProblem`, then solve with Nelder–Mead.
# The small initial simplex keeps the frequency search near each starting guess.

simplex = OptimizationOptimJL.Optim.AffineSimplexer(a=0.0, b=0.0005)
algorithm = NelderMead(initial_simplex=simplex)
after = Dict{Symbol,Any}()
histories = Dict{Symbol,Vector{Float64}}()

for (name, initial_label, target_label) in transfers
    initial = model.states[initial_label]
    target = model.states[target_label]
    history = Float64[]
    objective = function (m, gate, states)
        result = populations(m, gate, states...; samples=2)
        loss = 1 - last(result.target_population)
        push!(history, loss)
        return loss
    end

    println("Calibrating ", name)
    flush(stdout)
    problem = CalibrationProblem(model, trial_gates[name],
        ["drive/amplitude", "drive/carrier/frequency"], (initial, target); objective)
    solution = solve(problem, algorithm; maxiters=200, abstol=1e-5)
    gate = calibrated_gate(problem, solution.u)
    @assert SciMLBase.successful_retcode(solution.retcode)

    model.gates[name] = gate
    after[name] = populations(model, gate, initial, target)
    histories[name] = history
    @assert last(after[name].target_population) > 0.99
    @assert last(after[name].target_population) > last(before[name].target_population)
    println((; gate=name, amplitude_GHz=gate.parameters.drive.amplitude,
        frequency_GHz=gate.parameters.drive.carrier.frequency,
        before=last(before[name].target_population),
        after=last(after[name].target_population),
        final_leakage=last(after[name].leakage)))
end

# The reconstructed gate contains the calibrated parameters. We attach it to the model,
# then sample its full trajectory to check the result. The original trial gates
# remain available for comparison.
#
# ## 5. Compare the results
#
# Plot the trial and calibrated target populations for each transition:

results = Figure(size=(750, 850))
for (row, (name, initial, target)) in enumerate(transfers)
    ax = Axis(results[row, 1]; xlabel="Time (ns)", ylabel="Target population",
        title="$(initial) → $(target): $(name)")
    lines!(ax, before[name].times, before[name].target_population;
        label="Trial", color=:gray)
    lines!(ax, after[name].times, after[name].target_population;
        label="Calibrated", color=:purple)
    ylims!(ax, 0, 1.05)
    axislegend(ax; position=:lt)
end
results

# Each calibrated transfer should exceed 99%. The printed values also show the
# final population outside the chosen two-state pair. These are closed-system
# population transfers; phase errors and dissipation are not included in this
# calibration objective.
#
# We can also inspect how the best loss improved during calibration:

convergence = Figure(size=(700, 320))
ax_loss = Axis(convergence[1, 1]; xlabel="Objective evaluation",
    ylabel="Best transfer loss", yscale=log10)
for (name, _, _) in transfers
    history = histories[name]
    lines!(ax_loss, eachindex(history), max.(accumulate(min, history), 1e-12);
        label=string(name))
end
axislegend(ax_loss)
convergence

# ## 6. Save the result
#
# Save the calibrated model and the figures to a new temporary directory:

output_dir = mktempdir(; cleanup=false, prefix="transmon-resonator-")
model_path = QuantumDevices.save(joinpath(output_dir, "model"), model)
CairoMakie.save(joinpath(output_dir, "controls.png"), controls)
CairoMakie.save(joinpath(output_dir, "transfer.png"), results)
CairoMakie.save(joinpath(output_dir, "calibration.png"), convergence)
println("Saved model and figures to ", output_dir)
