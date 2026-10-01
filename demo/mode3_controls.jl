# Transmon–resonator control tutorial
# Run from the repository root:
#   julia --project=demo -e 'using Pkg; Pkg.instantiate()'
#   julia --project=demo demo/mode3_controls.jl [output_directory]
# Figures and a model/ bundle are saved in a new temporary directory by default.
# Include this file to use the helpers without running or saving the tutorial.

module Mode3ControlsDemo

using CairoMakie
using LinearAlgebra
using QuantumDevices
using QuantumToolbox
using SciMLBase

export TARGET, QUBIT_SETTINGS, SIDEBAND_SETTINGS, sideband_envelope, drive_pulse,
       make_gate, make_mode3_model, sideband_gap, sideband_population, plot_controls

# 1. Choose model and pulse parameters. Energies, frequencies, and amplitudes
# are in cycles/ns (GHz); time is in ns. Hamiltonians use angular units.
const TARGET = (
    EC=0.10283303447280807, EJ=26.96976142643705,
    resonator_frequency=6.228083962082612, g=0.026877206812551357,
    ng=0.0, n_cutoff=60, transmon_levels=10, resonator_levels=10,
)

const QUBIT_SETTINGS = (
    q_ge_0=(duration=92.96875, epsilon=0.00538, frequency=4.6040301557231,
        sigma=23.2421875, center=46.484375),
    q_ef_0=(duration=92.96875, epsilon=0.00385, frequency=4.49559986048746,
        sigma=23.2421875, center=46.484375),
)

const SIDEBAND_SETTINGS = (
    duration=205.05411739244443, epsilon=0.735,
    frequency=-2.869962347427361 + 0.037664663557229916,
    ramp_time=11.6257, k=2.0,
)

# 2. Define the envelope and carrier separately. GenericPulseFunction passes
# the owning gate duration to the callable as p.duration.
struct BumpSource end

function (::BumpSource)(p, t)
    x = 2t - 1
    abs(x) >= 1 && return 0.0
    x == 0 && return 1.0
    return exp(p.k * x^2 / (x^2 - 1))
end

"""Build a unit-height flattop using a smooth bump ramp."""
function sideband_envelope(settings=SIDEBAND_SETTINGS)
    bump = GenericPulseFunction(BumpSource(), (; k=settings.k))
    return ramped_flattop_pulse(1.0, settings.ramp_time; ramp=bump)
end

struct Carrier{E}
    envelope::E
end

function (carrier::Carrier)(p, t)
    return 2pi * p.epsilon * pulse_value(carrier.envelope, t, p.duration) *
           sinpi(2 * p.frequency * t)
end

"""Build `2π ε envelope(t) sin(2π frequency t)` in rad/ns."""
function drive_pulse(settings, envelope)
    return GenericPulseFunction(Carrier(envelope),
        (; epsilon=settings.epsilon, frequency=settings.frequency))
end

# 3. Attach a time-dependent coefficient to a symbolic charge operator.
# param(:drive) reads the matching entry in the gate parameters.
"""Attach a shaped charge drive to the transmon."""
function make_gate(settings, envelope)
    drive = drive_pulse(settings, envelope)
    return DeviceGate((; drive), param(:drive) * op(:q_charge), settings.duration)
end

# 4. Construct components, couple them, and select the retained dimensions.
# op(:q_charge) and op(:r_a) use the component names as prefixes.
"""
Build a transmon–resonator model and attach two Gaussian drives and a sideband.
The default 10×10 retained space resolves the sideband transfer. Smaller spaces
are useful for exploring model construction.
"""
function make_mode3_model(; params=TARGET, transmon_levels=params.transmon_levels,
                          resonator_levels=params.resonator_levels)
    q = make_transmon("q", 2pi * params.EC, 2pi * params.EJ,
        2params.n_cutoff + 1; ng=params.ng)
    r = make_resonator("r", 2pi * params.resonator_frequency, resonator_levels)
    interaction = param(:g) * op(:q_charge) *
                  (1im * (op(:r_a) - op(:r_adag)))
    model = make_model([q, r], interaction, (; g=2pi * params.g);
        truncation_dimensions=Dict(q => transmon_levels, r => resonator_levels),
        max_dimension=transmon_levels * resonator_levels)
    for (name, settings) in pairs(QUBIT_SETTINGS)
        envelope = gaussian_pulse(1.0, settings.sigma; center=settings.center)
        model.gates[name] = make_gate(settings, envelope)
    end
    model.gates[:sb_f0g1] = make_gate(SIDEBAND_SETTINGS, sideband_envelope())
    return model
end

# 5. Convert the model and gate into a numerical Hamiltonian, evolve a dressed
# state, and measure its overlap with the target state. Labels are zero-based.
"""Return the dressed `|f,0⟩`–`|g,1⟩` gap in cycles/ns."""
sideband_gap(model) = (model.others[(2, 0)] - model.others[(0, 1)]) / 2pi

"""Evolve `|f,0⟩` under the sideband and return times and `|g,1⟩` populations."""
function sideband_population(model; samples=101, abstol=1e-7, reltol=1e-7)
    gate = model.gates[:sb_f0g1]
    times = range(0, gate.duration; length=samples)
    result = sesolve(numerical(model, gate), model.states[(2, 0)], times;
        progress_bar=false, abstol, reltol)
    SciMLBase.successful_retcode(result.retcode) || error("Sideband solve did not succeed")
    target = model.states[(0, 1)]
    return times, [abs2(dot(target, state)) for state in result.states]
end

# 6. Inspect envelopes and carriers before interpreting the evolution.
"""Plot the qubit and f0g1 envelope and carrier coefficients."""
function plot_controls(model=make_mode3_model())
    fig = Figure(size=(1050, 760))
    q_env = Axis(fig[1, 1]; xlabel="Time (ns)", ylabel="Envelope amplitude (GHz)",
        title="Qubit Gaussian pulses")
    q_carrier = Axis(fig[1, 2]; xlabel="Time (ns)", ylabel="Drive coefficient / 2π (GHz)",
        title="Qubit carrier near pulse center")
    sb_env = Axis(fig[2, 1]; xlabel="Time (ns)", ylabel="Envelope amplitude (GHz)",
        title="f0g1 bump-ramp sideband")
    sb_carrier = Axis(fig[2, 2]; xlabel="Time (ns)", ylabel="Drive coefficient / 2π (GHz)",
        title="Sideband carrier near pulse center")

    for (name, color) in ((:q_ge_0, :royalblue), (:q_ef_0, :darkorange))
        gate = model.gates[name]
        T = gate.duration
        drive = gate.parameters.drive
        envelope = pulse_function(drive).envelope
        epsilon = drive.parameters.epsilon
        long = range(0, T; length=501)
        lines!(q_env, long, [epsilon * pulse_value(envelope, t, T) for t in long];
            label=string(name), color)
        short = range(T / 2 - 0.7, T / 2 + 0.7; length=501)
        lines!(q_carrier, short, [pulse_value(drive, t, T) / 2pi for t in short];
            label=string(name), color)
    end
    axislegend(q_env; position=:rt)
    axislegend(q_carrier; position=:rt)

    gate = model.gates[:sb_f0g1]
    T = gate.duration
    drive = gate.parameters.drive
    envelope = pulse_function(drive).envelope
    epsilon = drive.parameters.epsilon
    long = range(0, T; length=1001)
    short = range(T / 2 - 1.1, T / 2 + 1.1; length=701)
    lines!(sb_env, long, [epsilon * pulse_value(envelope, t, T) for t in long];
        color=:purple)
    lines!(sb_carrier, short, [pulse_value(drive, t, T) / 2pi for t in short];
        color=:purple)
    return fig
end

# 7. Save the instantiated model and its gates as a reusable bundle.
"""Run the tutorial and save its model bundle and figures to `output_dir`."""
function main(output_dir=mktempdir(; cleanup=false, prefix="mode3-tutorial-"))
    mkpath(output_dir)
    model = make_mode3_model()
    model_path = QuantumDevices.save(joinpath(output_dir, "model"), model)
    CairoMakie.save(joinpath(output_dir, "controls.png"), plot_controls(model))

    times, population = sideband_population(model; samples=201, abstol=1e-8, reltol=1e-8)
    fig = Figure(size=(700, 380))
    ax = Axis(fig[1, 1]; xlabel="Time (ns)", ylabel="|g,1⟩ population",
        title="f0g1 sideband transfer")
    lines!(ax, times, population; color=:purple, linewidth=3)
    ylims!(ax, 0, 1.05)
    CairoMakie.save(joinpath(output_dir, "transfer.png"), fig)
    println("Final |g,1⟩ population: ", last(population))
    println("Model bundle saved to ", model_path)
    println("Figures saved to ", abspath(output_dir))
    return (; model, model_path, times, population)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) ? Mode3ControlsDemo.main() : Mode3ControlsDemo.main(only(ARGS))
end
