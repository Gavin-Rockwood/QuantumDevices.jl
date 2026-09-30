# Transmon–resonator control tutorial
# Run from the repository root:
#   julia --project=demo -e 'using Pkg; Pkg.instantiate()'
#   julia --project=demo demo/mode3_controls.jl [figure_directory]
# With no directory argument, figures are saved to a new temporary directory.
# Include this file instead to work through its helpers interactively; inclusion
# defines the module without running the simulation or writing figures.

module Mode3ControlsDemo

using CairoMakie
using JSON3
using LinearAlgebra
using QuantumDevices
using QuantumToolbox
using SciMLBase

export load_controls, saved_envelope, drive_pulse, make_gate, make_mode3_model,
       sideband_gap, sideband_population, plot_controls

# 1. Load the saved control records. Frequencies and amplitudes are in GHz;
# times are in ns. Convert energy and drive coefficients to rad/ns with 2π.
const DATA_PATH = joinpath(@__DIR__, "data", "mode3_controls.json")

"""Read the relevant Mode3 records transcribed from the supplied save file."""
load_controls(path=DATA_PATH) = JSON3.read(read(path, String))

record(data, name) = data["Stuff"]["op_drive_params"][name]

# 2. Define the envelope and carrier separately. GenericPulseFunction passes
# the owning gate duration to the callable as p.duration.

# This normalized source is the old Bump_Envelope on the unit interval.
struct BumpSource end

function (::BumpSource)(p, t)
    x = 2t - 1
    abs(x) >= 1 && return 0.0
    x == 0 && return 1.0
    return exp(p.k * x^2 / (x^2 - 1))
end

"""Rebuild the saved Gaussian or Bump_Ramp envelope at unit peak height."""
function saved_envelope(pulse_record)
    shape = String(pulse_record["Envelope"])
    args = pulse_record["Envelope Args"]
    if shape == "Guassian" # Spelling retained from the source save format.
        return gaussian_pulse(1.0, Float64(args["sigma"]); center=Float64(args["mu"]))
    elseif shape == "Bump_Ramp"
        T = Float64(pulse_record["pulse_time"])
        isapprox(Float64(args["pulse_time"]), T; rtol=0, atol=1e-10) ||
            throw(ArgumentError("Bump_Ramp duration differs from pulse_time"))
        bump = GenericPulseFunction(BumpSource(), (; k=Float64(args["k"])))
        return ramped_flattop_pulse(1.0, Float64(args["ramp_time"]); ramp=bump)
    end
    throw(ArgumentError("Unsupported saved envelope: $shape"))
end

struct Carrier{E}
    envelope::E
end

function (carrier::Carrier)(p, t)
    return 2pi * p.epsilon * pulse_value(carrier.envelope, t, p.duration) *
           sinpi(2 * p.frequency * t)
end

"""
Rebuild `2π ε envelope(t) sin(2π(freq_d + shift)t)` in rad/ns.
The source save uses GHz and ns, and `shift` is a frequency offset in GHz.
"""
function drive_pulse(pulse_record)
    envelope = saved_envelope(pulse_record)
    epsilon = Float64(pulse_record["epsilon"])
    frequency = Float64(pulse_record["freq_d"] + pulse_record["shift"])
    return GenericPulseFunction(Carrier(envelope), (; epsilon, frequency))
end

# 3. Attach a time-dependent coefficient to a symbolic charge operator.
# param(:drive) reads the matching entry in the gate parameters.

"""Attach a saved control to the transmon charge operator."""
function make_gate(data, name)
    pulse_record = record(data, name)
    duration = Float64(pulse_record["pulse_time"])
    return DeviceGate((; drive=drive_pulse(pulse_record)),
        param(:drive) * op(:q_charge), duration)
end

# 4. Construct components, couple them, and select the retained dimensions.
# op(:q_charge) and op(:r_a) use the component names as prefixes.

"""
Build the saved transmon–Mode3 model. The default 10×10 retained space matches
the save; smaller `transmon_levels` and `resonator_levels` speed up a tutorial run.
The parent transmon charge basis always uses the saved `Nt_cut`.
"""
function make_mode3_model(data=load_controls(); transmon_levels=nothing,
                          resonator_levels=nothing)
    config = data["Main_Config"]
    length(config["E_oscs"]) == length(config["gs"]) == length(config["Nrs"]) == 1 ||
        throw(ArgumentError("This demo expects one resonator mode"))
    nt = transmon_levels === nothing ? Int(config["Nt"]) : Int(transmon_levels)
    nr = resonator_levels === nothing ? Int(config["Nrs"][1]) : Int(resonator_levels)
    q = make_transmon("q", 2pi * config["E_C"], 2pi * config["E_J"],
        2 * Int(config["Nt_cut"]) + 1; ng=config["ng"])
    r = make_resonator("r", 2pi * config["E_oscs"][1], Int(config["Nrs"][1]))
    # The source constructor uses charge ⊗ i(a-a†), not charge ⊗ (a+a†).
    interaction = param(:g) * op(:q_charge) *
                  (1im * (op(:r_a) - op(:r_adag)))
    model = make_model([q, r], interaction, (; g=2pi * config["gs"][1]);
        truncation_dimensions=Dict(q => nt, r => nr), max_dimension=nt * nr)
    for name in keys(data["Stuff"]["op_drive_params"])
        model.gates[Symbol(name)] = make_gate(data, name)
    end
    return model
end

# 5. Convert the model and gate into a numerical Hamiltonian, evolve a dressed
# state, and measure its overlap with the target state. Labels are zero-based.

"""Return the dressed `|f,0⟩`–`|g,1⟩` gap in cycles/ns."""
sideband_gap(model) = (model.others[(2, 0)] - model.others[(0, 1)]) / 2pi

"""Evolve `|f,0⟩` under the saved sideband and return times and `|g,1⟩` populations."""
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

"""Plot the saved qubit and f0g1 envelope and carrier coefficients."""
function plot_controls(data=load_controls())
    fig = Figure(size=(1050, 760))
    q_env = Axis(fig[1, 1]; xlabel="Time (ns)", ylabel="Envelope amplitude (GHz)",
        title="Qubit Gaussian pulses")
    q_carrier = Axis(fig[1, 2]; xlabel="Time (ns)", ylabel="Drive coefficient / 2π (GHz)",
        title="Qubit carrier near pulse center")
    sb_env = Axis(fig[2, 1]; xlabel="Time (ns)", ylabel="Envelope amplitude (GHz)",
        title="f0g1 bump-ramp sideband")
    sb_carrier = Axis(fig[2, 2]; xlabel="Time (ns)", ylabel="Drive coefficient / 2π (GHz)",
        title="Sideband carrier near pulse center")

    for (name, color) in (("q_ge_0", :royalblue), ("q_ef_0", :darkorange))
        saved = record(data, name)
        T = Float64(saved["pulse_time"])
        epsilon = Float64(saved["epsilon"])
        envelope = saved_envelope(saved)
        drive = drive_pulse(saved)
        lines!(q_env, range(0, T; length=501),
            [epsilon * pulse_value(envelope, t, T) for t in range(0, T; length=501)];
            label=name, color)
        short = range(T / 2 - 0.7, T / 2 + 0.7; length=501)
        lines!(q_carrier, short, [pulse_value(drive, t, T) / 2pi for t in short];
            label=name, color)
    end
    axislegend(q_env; position=:rt)
    axislegend(q_carrier; position=:rt)

    saved = record(data, "sb_f0g1")
    T = Float64(saved["pulse_time"])
    epsilon = Float64(saved["epsilon"])
    envelope = saved_envelope(saved)
    drive = drive_pulse(saved)
    long = range(0, T; length=1001)
    short = range(T / 2 - 1.1, T / 2 + 1.1; length=701)
    lines!(sb_env, long, [epsilon * pulse_value(envelope, t, T) for t in long];
        color=:purple)
    lines!(sb_carrier, short, [pulse_value(drive, t, T) / 2pi for t in short];
        color=:purple)
    return fig
end

"""Run the tutorial and save its two figures to `output_dir`."""
function main(output_dir=mktempdir(; cleanup=false, prefix="mode3-tutorial-"))
    mkpath(output_dir)
    controls = load_controls()
    CairoMakie.save(joinpath(output_dir, "controls.png"), plot_controls(controls))

    # Use the saved 10×10 space to reproduce the sideband transfer.
    model = make_mode3_model(controls)
    times, population = sideband_population(model; samples=201, abstol=1e-8, reltol=1e-8)
    fig = Figure(size=(700, 380))
    ax = Axis(fig[1, 1]; xlabel="Time (ns)", ylabel="|g,1⟩ population",
        title="f0g1 sideband transfer")
    lines!(ax, times, population; color=:purple, linewidth=3)
    ylims!(ax, 0, 1.05)
    CairoMakie.save(joinpath(output_dir, "transfer.png"), fig)
    println("Final |g,1⟩ population: ", last(population))
    println("Figures saved to ", abspath(output_dir))
    return (; model, times, population)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) ? Mode3ControlsDemo.main() : Mode3ControlsDemo.main(only(ARGS))
end
