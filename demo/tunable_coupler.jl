# Tunable coupler control tutorial
# Run from the repository root:
#   julia --project=demo -e 'using Pkg; Pkg.instantiate()'
#   julia --project=demo demo/tunable_coupler.jl [output_directory]
# Figures and a model/ bundle are saved in a new temporary directory by default.
# Include this file instead to work through its helpers interactively; inclusion
# defines the module without running the simulation or writing figures.
# For a step-by-step calibration and projected gate comparison, see
# docs/src/tutorials/tunable_coupler_control.jl.

module TunableCouplerDemo

using CairoMakie
using LinearAlgebra
using QuantumDevices
using QuantumToolbox
using SciMLBase

export TARGET, SWAP_SETTINGS, QUARTER_SETTINGS, make_coupler_model,
       swap_gate, quarter_gate, plot_controls, swap_population

# 1. Choose the circuit and pulse parameters. Flux amplitudes use model phi.
# Energies and drive frequencies are in cycles/ns; durations are in ns.
const TARGET = (
    EC1=0.198, EC2=0.18, ECc=0.097,
    EJ1=11.143, EJ2=11.15, EJc=33.528,
    EC1c=4.3, EC2c=4.5, EC12=161.4,
)

const SWAP_SETTINGS = (
    duration=28.135660302912168,
    ramp=1.8727503955158111,
    q1_flux=0.1318136929161456,
    coupler_flux=0.3873022470423121,
)

const QUARTER_SETTINGS = (
    q1=(frequency=3.98687197160622, epsilon=0.024775449345990953,
        phase=1.5501484428191823, duration=9.027992059968785),
    q2=(frequency=3.813857655541284, epsilon=0.019776450971418814,
        phase=-4.73026239790409, duration=11.010965386740736),
)

# 2. Construct three transmons and couple their charge operators.
# Component names set the prefixes used in symbolic operators and parameters.

"""Build the three-mode transmon–transmon–coupler model with capacitive charge coupling."""
function make_coupler_model(; params=TARGET, retained_levels=3, n_cutoff=60)
    n_cutoff isa Integer && n_cutoff > 0 || throw(ArgumentError("n_cutoff must be positive"))
    retained_levels isa Integer && 1 <= retained_levels <= 2n_cutoff+1 ||
        throw(ArgumentError("retained_levels must fit the parent charge basis"))
    # Split the total Josephson energy equally between the two junctions.
    q1 = make_tunable_transmon("q1", params.EC1, params.EJ1/2,
        params.EJ1/2, 2n_cutoff+1)
    q2 = make_tunable_transmon("q2", params.EC2, params.EJ2/2,
        params.EJ2/2, 2n_cutoff+1)
    c = make_tunable_transmon("c", params.ECc, params.EJc/2,
        params.EJc/2, 2n_cutoff+1)

    η = params.EC12 * params.ECc / (params.EC1c * params.EC2c)
    g1c = 8 * params.EC1 * params.ECc / params.EC1c
    g2c = 8 * params.EC2 * params.ECc / params.EC2c
    g12 = 8 * (1 + η) * params.EC1 * params.EC2 / params.EC12
    interaction = param(:g1c)*op(:q1_charge)*op(:c_charge) +
                  param(:g2c)*op(:q2_charge)*op(:c_charge) +
                  param(:g12)*op(:q1_charge)*op(:q2_charge)
    model = make_model([q1, q2, c], interaction,
        (; g1c, g2c, g12);
        truncation_dimensions=Dict(q1=>retained_levels, q2=>retained_levels,
                                   c=>retained_levels),
        max_dimension=retained_levels^3)
    model.gates[:swap] = swap_gate()
    for axis in (:X, :Y), qubit in (1, 2)
        model.gates[Symbol("quarter_", lowercase(string(axis)), qubit)] =
            quarter_gate(qubit, axis)
    end
    return model
end

# 3. Change existing model parameters with a flux gate. A zero additional
# Hamiltonian is sufficient: q1_phi and c_phi already occur in the model.

"""Exchange pulse with sin² flux ramps on q1 and the coupler."""
function swap_gate(; settings=SWAP_SETTINGS)
    T, ramp = settings.duration, settings.ramp
    q1_flux = Pulse(RampedFlattop(ramp); amplitude=settings.q1_flux, duration=T)
    # The coupler starts one ramp late and ends one ramp early.
    c_flux = Pulse(RampedFlattop(ramp); amplitude=settings.coupler_flux,
        duration=T-2ramp, delay=ramp)
    return DeviceGate((; q1_phi=q1_flux, c_phi=c_flux), 0)
end

# 4. Add microwave drives with reusable envelopes and carriers.

"""Build one of the quarter-X/Y charge drives on qubit 1 or 2."""
function quarter_gate(qubit::Integer, axis::Symbol)
    qubit in (1, 2) || throw(ArgumentError("qubit must be 1 or 2"))
    axis in (:X, :Y) || throw(ArgumentError("axis must be :X or :Y"))
    settings = qubit == 1 ? QUARTER_SETTINGS.q1 : QUARTER_SETTINGS.q2
    phase = axis === :Y ? settings.phase - 3pi/2 : settings.phase
    drive = Pulse(SineSquared(); amplitude=settings.epsilon, duration=settings.duration,
        carrier=SineCarrier(settings.frequency; phase))
    return DeviceGate((; drive), param(:drive)*op(Symbol("q", qubit, "_charge")))
end

# 5. Evolve a dressed state and measure source and target populations.
# State labels follow the component order [q1, q2, c] and are zero-based.

"""Return the dressed |100⟩ and |010⟩ populations under the flux pulse."""
function swap_population(model; samples=121, abstol=1e-8, reltol=1e-8)
    gate = model.gates[:swap]
    times = range(0, gate.duration; length=samples)
    result = sesolve(model, gate, model.states[(1, 0, 0)], times;
        progress_bar=false, abstol, reltol)
    SciMLBase.successful_retcode(result.retcode) || error("Flux-pulse solve failed")
    source = model.states[(1, 0, 0)]
    target = model.states[(0, 1, 0)]
    amplitudes = state_amplitudes(Dict(:source => source, :target => target), result)
    return times, abs2.(amplitudes[:source]), abs2.(amplitudes[:target])
end

# 6. Inspect the flux excursions and microwave coefficients.

"""Plot the flux excursions and quarter-X/Y drive coefficients."""
function plot_controls(model=make_coupler_model())
    fig = Figure(size=(980, 690))
    flux_axis = Axis(fig[1, 1:2]; xlabel="Time (ns)", ylabel="Flux (Φ₀)",
        title="Coupler exchange pulse")
    T = model.gates[:swap].duration
    ts = range(0, T; length=501)
    gate = model.gates[:swap]
    lines!(flux_axis, ts, [gate.parameters.q1_phi(t) for t in ts];
        color=:royalblue, label="Qubit 1")
    lines!(flux_axis, ts, [gate.parameters.c_phi(t) for t in ts];
        color=:purple, label="Coupler")
    axislegend(flux_axis; position=:rt)
    for qubit in (1, 2)
        ax = Axis(fig[2, qubit]; xlabel="Time (ns)",
            ylabel="Charge-drive coefficient (GHz)",
            title="Qubit $qubit quarter rotations")
        for (axis, color) in ((:X, :royalblue), (:Y, :darkorange))
            gate = model.gates[Symbol("quarter_", lowercase(string(axis)), qubit)]
            T = gate.duration
            ts = range(0, T; length=1001)
            lines!(ax, ts, [gate.parameters.drive(t) for t in ts];
                color, label="quarter-$axis")
        end
        axislegend(ax; position=:rt)
    end
    return fig
end

# 7. Save the instantiated model and its gates as a reusable bundle.
"""Run the tutorial and save its model bundle and figures to `output_dir`."""
function main(output_dir=mktempdir(; cleanup=false, prefix="coupler-tutorial-"))
    mkpath(output_dir)
    model = make_coupler_model()
    model_path = QuantumDevices.save(joinpath(output_dir, "model"), model)
    CairoMakie.save(joinpath(output_dir, "controls.png"), plot_controls(model))

    times, p100, p010 = swap_population(model)
    fig = Figure(size=(700, 370))
    ax = Axis(fig[1, 1]; xlabel="Time (ns)", ylabel="Dressed-state population",
        title="Coupler flux pulse")
    lines!(ax, times, p100; label="|100⟩", color=:royalblue, linewidth=2.5)
    lines!(ax, times, p010; label="|010⟩", color=:purple, linewidth=2.5)
    axislegend(ax; position=:rc)
    ylims!(ax, 0, 1.05)
    CairoMakie.save(joinpath(output_dir, "transfer.png"), fig)
    println("Final |010⟩ population: ", last(p010))
    println("Model bundle saved to ", model_path)
    println("Figures saved to ", abspath(output_dir))
    return (; model, model_path, times, p100, p010)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) ? TunableCouplerDemo.main() : TunableCouplerDemo.main(only(ARGS))
end
