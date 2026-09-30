# Tunable coupler control tutorial
# Run from the repository root:
#   julia --project=demo -e 'using Pkg; Pkg.instantiate()'
#   julia --project=demo demo/tunable_coupler.jl [figure_directory]
# With no directory argument, figures are saved to a new temporary directory.
# Include this file instead to work through its helpers interactively; inclusion
# defines the module without running the simulation or writing figures.

module TunableCouplerDemo

using CairoMakie
using LinearAlgebra
using QuantumDevices
using QuantumToolbox
using SciMLBase

export TARGET, SWAP_SETTINGS, QUARTER_SETTINGS, make_coupler_model,
       swap_gate, quarter_gate, plot_controls, swap_population

# 1. Choose the circuit and pulse parameters.
# Energies and drive frequencies are in cycles/ns; durations are in ns.
const TARGET = (
    EC1=0.198, EC2=0.18, ECc=0.097,
    EJ1=11.143, EJ2=11.15, EJc=33.528,
    EC1c=4.3, EC2c=4.5, EC12=161.4,
)

const SWAP_SETTINGS = (
    duration=28.135660302912168,
    ramp=1.8727503955158111,
    q1_flux=0.06589740537034883,
    coupler_flux=0.1941640246800509,
)

const QUARTER_SETTINGS = (
    q1=(frequency=3.98687197160622, epsilon=0.024775449345990953,
        phase=1.5501484428191823, duration=9.027992059968785),
    q2=(frequency=3.813857655541284, epsilon=0.019776450971418814,
        phase=-4.73026239790409, duration=11.010965386740736),
)

# 2. Construct three transmons and couple their charge operators.
# Component names set the prefixes used in symbolic operators and parameters.

"""Build the three-mode transmon–transmon–coupler model from the supplied values."""
function make_coupler_model(; params=TARGET, retained_levels=3, n_cutoff=60)
    n_cutoff isa Integer && n_cutoff > 0 || throw(ArgumentError("n_cutoff must be positive"))
    retained_levels isa Integer && 1 <= retained_levels <= 2n_cutoff+1 ||
        throw(ArgumentError("retained_levels must fit the parent charge basis"))
    # Equal junction energies give d=0. Package phi=2Φ reproduces cos(2πΦ).
    q1 = make_tunable_transmon("q1", 2pi*params.EC1, pi*params.EJ1,
        pi*params.EJ1, 2n_cutoff+1)
    q2 = make_tunable_transmon("q2", 2pi*params.EC2, pi*params.EJ2,
        pi*params.EJ2, 2n_cutoff+1)
    c = make_tunable_transmon("c", 2pi*params.ECc, pi*params.EJc,
        pi*params.EJc, 2n_cutoff+1)

    η = params.EC12 * params.ECc / (params.EC1c * params.EC2c)
    g1c = 8 * params.EC1 * params.ECc / params.EC1c
    g2c = 8 * params.EC2 * params.ECc / params.EC2c
    g12 = 8 * (1 + η) * params.EC1 * params.EC2 / params.EC12
    interaction = param(:g1c)*op(:q1_charge)*op(:c_charge) +
                  param(:g2c)*op(:q2_charge)*op(:c_charge) +
                  param(:g12)*op(:q1_charge)*op(:q2_charge)
    model = make_model([q1, q2, c], interaction,
        (; g1c=2pi*g1c, g2c=2pi*g2c, g12=2pi*g12);
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

"""Flux pulse labeled `SWAP` in the source, with linear ramps on q1 and coupler."""
function swap_gate(; settings=SWAP_SETTINGS)
    T, ramp = settings.duration, settings.ramp
    q1_flux = ramped_flattop_pulse(2 * settings.q1_flux, ramp;
        ramp=:linear, stop=T)
    # The coupler pulse begins one ramp late and ends one ramp early.
    c_flux = ramped_flattop_pulse(2 * settings.coupler_flux, ramp;
        ramp=:linear, start=ramp, stop=T-ramp)
    return DeviceGate((; q1_phi=q1_flux, c_phi=c_flux), 0, T)
end

# 4. Add microwave drives with a separate envelope and carrier. The gate
# supplies p.duration; GenericPulseFunction stores only the other parameters.

struct ChargeCarrier{E}
    envelope::E
end

function (carrier::ChargeCarrier)(p, t)
    return 2pi*p.epsilon * pulse_value(carrier.envelope, t, p.duration) *
           sin(2pi*p.frequency*t + p.phase)
end

"""Build one of the saved quarter-X/Y charge drives on qubit 1 or 2."""
function quarter_gate(qubit::Integer, axis::Symbol)
    qubit in (1, 2) || throw(ArgumentError("qubit must be 1 or 2"))
    axis in (:X, :Y) || throw(ArgumentError("axis must be :X or :Y"))
    settings = qubit == 1 ? QUARTER_SETTINGS.q1 : QUARTER_SETTINGS.q2
    phase = axis === :Y ? settings.phase - 3pi/2 : settings.phase
    envelope = sine_squared_pulse(1.0, 1/(2settings.duration))
    drive = GenericPulseFunction(ChargeCarrier(envelope),
        (; epsilon=settings.epsilon, frequency=settings.frequency, phase))
    return DeviceGate((; drive), param(:drive)*op(Symbol("q", qubit, "_charge")),
        settings.duration)
end

# 5. Evolve a dressed state and measure source and target populations.
# State labels follow the component order [q1, q2, c] and are zero-based.

"""Return the dressed |100⟩ and |010⟩ populations under the saved flux pulse."""
function swap_population(model; samples=121, abstol=1e-8, reltol=1e-8)
    gate = model.gates[:swap]
    times = range(0, gate.duration; length=samples)
    result = sesolve(numerical(model, gate), model.states[(1, 0, 0)], times;
        progress_bar=false, abstol, reltol)
    SciMLBase.successful_retcode(result.retcode) || error("Flux-pulse solve failed")
    source = model.states[(1, 0, 0)]
    target = model.states[(0, 1, 0)]
    return times,
        [abs2(dot(source, state)) for state in result.states],
        [abs2(dot(target, state)) for state in result.states]
end

# 6. Inspect the flux excursions and microwave coefficients.

"""Plot the saved flux excursions and quarter-X/Y drive coefficients."""
function plot_controls(model=make_coupler_model())
    fig = Figure(size=(980, 690))
    flux_axis = Axis(fig[1, 1:2]; xlabel="Time (ns)", ylabel="Flux (Φ₀)",
        title="Coupler exchange pulse")
    T = model.gates[:swap].duration
    ts = range(0, T; length=501)
    gate = model.gates[:swap]
    lines!(flux_axis, ts, [pulse_value(gate.parameters.q1_phi, t, T)/2 for t in ts];
        color=:royalblue, label="Qubit 1")
    lines!(flux_axis, ts, [pulse_value(gate.parameters.c_phi, t, T)/2 for t in ts];
        color=:purple, label="Coupler")
    axislegend(flux_axis; position=:rt)
    for qubit in (1, 2)
        ax = Axis(fig[2, qubit]; xlabel="Time (ns)",
            ylabel="Charge-drive coefficient / 2π (GHz)",
            title="Qubit $qubit quarter rotations")
        for (axis, color) in ((:X, :royalblue), (:Y, :darkorange))
            gate = model.gates[Symbol("quarter_", lowercase(string(axis)), qubit)]
            T = gate.duration
            ts = range(0, T; length=1001)
            lines!(ax, ts, [pulse_value(gate.parameters.drive, t, T)/2pi for t in ts];
                color, label="quarter-$axis")
        end
        axislegend(ax; position=:rt)
    end
    return fig
end

"""Run the tutorial and save its two figures to `output_dir`."""
function main(output_dir=mktempdir(; cleanup=false, prefix="coupler-tutorial-"))
    mkpath(output_dir)
    model = make_coupler_model()
    CairoMakie.save(joinpath(output_dir, "controls.png"), plot_controls(model))

    times, p100, p010 = swap_population(model)
    fig = Figure(size=(700, 370))
    ax = Axis(fig[1, 1]; xlabel="Time (ns)", ylabel="Dressed-state population",
        title="Coupler flux pulse")
    lines!(ax, times, p100; label="|100⟩", color=:royalblue, linewidth=2.5)
    lines!(ax, times, p010; label="|010⟩", color=:purple, linewidth=2.5)
    axislegend(ax; position=:rt)
    ylims!(ax, 0, 1.05)
    CairoMakie.save(joinpath(output_dir, "transfer.png"), fig)
    println("Final |010⟩ population: ", last(p010))
    println("Figures saved to ", abspath(output_dir))
    return (; model, times, p100, p010)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) ? TunableCouplerDemo.main() : TunableCouplerDemo.main(only(ARGS))
end
