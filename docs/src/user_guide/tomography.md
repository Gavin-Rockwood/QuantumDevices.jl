# Model tomography

First fit `EC`, `EJ`, the coupling `g`, and cavity frequency from spectroscopy.
Then fit an extra `chi*N` term from gate-sequence measurements after a TLS flip.
Frequencies are in GHz and times in ns. The example uses synthetic data; replace
the observations with your measurements to fit a device.

This example is self-contained. Install OptimizationOptimJL in your application
environment to run the fits. Start Julia with `--threads=4` to build distinct
gate propagators in parallel; a single-thread session also works.

```@example cavity_tomography
using QuantumDevices, QuantumToolbox, OptimizationOptimJL, SciMLBase, CairoMakie

CairoMakie.activate!(type="svg")
threaded = true; # Use a kernel started with --threads=4 for parallel gate builds.
```

## 1. Fit the baseline model

The four inputs are the transmon g-e and e-f frequencies, the cavity 0→1
frequency with the transmon in g, and the conditional cavity shift
`(e1 − e0) − (g1 − g0)`. Here `(j, n)` labels a dressed state tracked from
transmon level `j` and cavity occupation `n`. Transmon anharmonicity
`f_ef − f_ge` is a derived quantity, rather than another independent input.

`ModelDatum` compares a scalar prediction against an observed value. Its residual
is `(prediction - value) / uncertainty`; `uncertainty=0.001` below means 1 MHz.
The scalar fit uses the sum of squared residuals and does not evolve any gates.
`real(...)` removes imaginary numerical roundoff from the Hermitian spectrum.

```@example cavity_tomography
function cavity_transmon_model(EC, EJ, g, cavity_frequency)
    transmon = make_transmon("t", EC, EJ, 17)
    cavity = make_resonator("c", cavity_frequency, 3)
    interaction = param(:g)*op(:t_charge)*(op(:c_a) + op(:c_adag))
    make_model([transmon, cavity], interaction, (; g);
        truncation_dimensions=Dict(transmon => 4))
end

quantities = [
    m -> real(m.spectrum[(1, 0)] - m.spectrum[(0, 0)]),
    m -> real(m.spectrum[(2, 0)] - m.spectrum[(1, 0)]),
    m -> real(m.spectrum[(0, 1)] - m.spectrum[(0, 0)]),
    m -> real((m.spectrum[(1, 1)] - m.spectrum[(1, 0)]) -
              (m.spectrum[(0, 1)] - m.spectrum[(0, 0)])),
]

# Replace these four values with measured frequencies, in GHz.
reference = cavity_transmon_model(0.2, 10.0, 0.08, 4.5)
measured = (
    ge=quantities[1](reference),
    ef=quantities[2](reference),
    cavity_01=quantities[3](reference),
    conditional_cavity_shift=quantities[4](reference),
)
spectral_data = [ModelDatum(f, value; uncertainty=0.001)
    for (f, value) in zip(quantities, values(measured))];
```

```@example cavity_tomography
initial = cavity_transmon_model(0.21, 9.7, 0.07, 4.45)
paths = ["components/1/EC", "components/1/EJ",
         "coupling_parameters/g", "components/2/frequency"]
spectral_problem = ModelTomographyProblem(initial, spectral_data, paths)
spectral_solution = solve(spectral_problem, OptimizationOptimJL.NelderMead();
    maxiters=1500, g_abstol=1e-10)
@assert SciMLBase.successful_retcode(spectral_solution.retcode)
baseline = reconstructed_model(spectral_problem, spectral_solution.u)
for (name, value) in zip(("EC", "EJ", "g", "Cavity frequency"), spectral_solution.u)
    println(name, " = ", round(value; digits=6), " GHz")
end
println("Transmon anharmonicity = ", round(1000*(measured.ef - measured.ge); digits=3), " MHz")
@assert isapprox(spectral_solution.u, [0.2, 10.0, 0.08, 4.5]; rtol=1e-3) # hide
```

## 2. Add the TLS-induced shift

Keep the fitted physical parameters fixed and add `chi*N`, where `N` counts
transmon energy-level excitations. This models a constant offset after a TLS
flip, rather than a TLS switching during a gate sequence.

The following controls are laboratory-frame microwave pulses on the transmon
charge operator. The readout pulse phase follows the same oscillator clock
through the wait and includes a quadrature offset to distinguish the shift sign.

```@example cavity_tomography
transmon = baseline.components[1]
transmon_with_N = Component(transmon.name, transmon.parameters,
    merge(transmon.operators, (; N=num(transmon.dimension))),
    transmon.hamiltonian, transmon.type, transmon.dimension)
shift_trial = make_model([transmon_with_N, baseline.components[2]],
    baseline.interactions + param(:chi)*op(:t_N),
    merge(baseline.coupling_parameters, (; chi=0.005));
    truncation_dimensions=Dict(transmon_with_N => 4))

carrier = quantities[1](baseline)
pulse_time = 5.0
idle_time = 4.0
wait_counts = 0:6
prepare_pulse = Pulse(Constant(); duration=pulse_time, amplitude=0.05,
    carrier=SineCarrier(carrier; phase=pi/2))
shift_trial.gates["prepare"] = DeviceGate((drive=prepare_pulse,), param(:drive)*op(:t_charge), pulse_time)
shift_trial.gates["idle"] = DeviceGate((;), val(0), idle_time)
for k in wait_counts
    start_time = pulse_time + k*idle_time
    readout_pulse = Pulse(Constant(); duration=pulse_time, amplitude=0.05,
        carrier=SineCarrier(carrier; phase=pi + 2pi*carrier*start_time))
    shift_trial.gates["readout_$k"] = DeviceGate((drive=readout_pulse,), param(:drive)*op(:t_charge), pulse_time)
end
```

## 3. Fit the shift from sequence data

Prepare the dressed ground state, apply a pulse, wait, and apply the readout
pulse. Measure all four retained transmon levels, tracing out the cavity.
All outcomes enter the fit, including leakage above the g-e subspace.

Each `SequenceData` record contains its gate sequence, preparation, measurement,
and observed distribution. The preparation function is evaluated against each
candidate model. `Measurement(num(4); acts_on=[1])` orders outcomes as
transmon levels 0, 1, 2, 3, so index 2 below is the level-1 probability.

```@example cavity_tomography
true_chi = 0.015  # Synthetic post-flip shift: 15 MHz.
shift_target = setpath(shift_trial, "coupling_parameters/chi", true_chi)
preparation = m -> m.states[(0, 0)]
measurement = Measurement(num(4); acts_on=[1])
sequences = [vcat(["prepare"], fill("idle", k), ["readout_$k"]) for k in wait_counts]
# Templates allow all synthetic predictions to share each gate propagator.
records = [SequenceData(s, preparation, measurement; probabilities=[1.0, 0, 0, 0])
    for s in sequences]
observed = predict_sequences(shift_target, records; threaded)
# Replace observed with your post-flip measurement distributions.
sequence_data = [SequenceData(s, preparation, measurement; probabilities=p)
    for (s, p) in zip(sequences, observed)]

shift_problem = ModelTomographyProblem(shift_trial, sequence_data,
    ["coupling_parameters/chi"]; threaded)
shift_solution = solve(shift_problem, OptimizationOptimJL.NelderMead();
    maxiters=200, g_abstol=1e-12)
@assert SciMLBase.successful_retcode(shift_solution.retcode)
fitted = reconstructed_model(shift_problem, shift_solution.u)
println("True shift:   ", 1000*true_chi, " MHz")
println("Fitted shift: ", round(1000*shift_solution.u[1]; digits=4), " MHz")
@assert isapprox(shift_solution.u[1], true_chi; atol=1e-5) # hide
```

## Compare the predictions

The plot shows level-1 population. The fit uses the complete distributions;
it does not condition on remaining within the first two transmon levels.

```@example cavity_tomography
before_fit = predict_sequences(shift_trial, sequence_data; threaded)
after_fit = predict_sequences(fitted, sequence_data; threaded)
wait_times = idle_time .* collect(wait_counts)
figure = Figure(size=(700, 400))
axis = Axis(figure[1, 1]; xlabel="Waiting time (ns)",
    ylabel="Transmon level-1 probability", title="Fit a post-flip frequency shift")
scatter!(axis, wait_times, [p[2] for p in observed];
    label="Observed", color=:black, markersize=12)
lines!(axis, wait_times, [p[2] for p in before_fit];
    label="Initial shift", color=:darkorange, linestyle=:dash)
lines!(axis, wait_times, [p[2] for p in after_fit];
    label="Fitted shift", color=:dodgerblue)
ylims!(axis, 0, 1)
axislegend(axis; position=:rt)
figure
```

## Using experimental data

For experimental data, use measured uncertainties in the spectral fit and
`counts` (or probabilities with `shots`) in the sequence fit. Check other initial
guesses and additional sequences before accepting a result. This example starts
near a positive `g`; spectroscopy may determine only its magnitude. Check basis
convergence and dressed-state labels near avoided crossings. The added `N`
operator is defined after the physical fit, so `EC` and `EJ` stay fixed in stage 2.

`ModelTomographyProblem` accepts scalar observations or sequence records, not a
mixture. Pass the reconstructed spectral model into the sequence stage as shown
here. Every candidate is rebuilt from its original trial model; fitting does not
overwrite that model. Gate-only updates reuse idle physics, interaction updates
retrack dressed states, and component updates such as `EC` rebuild the affected
component and dependent projections.

The evolution here is unitary: decoherence and switching during acquisition are
not included. See the [API reference](../reference/tomography.md) for measurement
options, loss functions, and solver configuration.
