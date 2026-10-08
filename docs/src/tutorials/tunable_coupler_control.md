```@meta
EditURL = "tunable_coupler_control.jl"
```

# Tunable coupler control

This tutorial constructs two transmons coupled through a third, flux-tunable
transmon. A flux pulse brings the qubits into resonance and changes their
effective coupling, allowing an excitation to move between them:

$$|1,0,0\rangle \longleftrightarrow |0,1,0\rangle.$$

We define the pulse, simulate the transfer, and calibrate its flux amplitudes.
We then project the evolution onto the computational states to inspect phases
and leakage as well as populations.

## 1. Set up the environment

Run the code blocks in order, or run this source file from the repository root:

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/src/tutorials/tunable_coupler_control.jl
```

Time is in **ns**, energies and Hamiltonian coefficients are in **GHz**, flux
is in units of the flux quantum **Φ₀**, and carrier phases are in radians.
The model-aware `sesolve` and projected `numerical` methods apply the `2π`
conversion for evolution internally.

```@example coupler
using QuantumDevices, QuantumToolbox
using OptimizationOptimJL
```

## 2. Build the coupled model

Each component is a symmetric split-junction transmon. Its effective Josephson
energy is $E_J(\phi)=E_{J,\max}|\cos(\pi\phi)|$. Split the maximum Josephson
energy equally between the two junctions when calling `make_tunable_transmon`.
The charge operators couple the two qubits to the coupler and to each other.

```@example coupler
EC1, EC2, ECc = 0.198, 0.18, 0.097
EJ1, EJ2, EJc = 11.143, 11.15, 33.528
EC1c, EC2c, EC12 = 4.3, 4.5, 161.4
charge_cutoff = 60
retained_levels = 3

q1 = make_tunable_transmon("q1", EC1, EJ1 / 2, EJ1 / 2, 2charge_cutoff + 1)
q2 = make_tunable_transmon("q2", EC2, EJ2 / 2, EJ2 / 2, 2charge_cutoff + 1)
coupler = make_tunable_transmon("c", ECc, EJc / 2, EJc / 2, 2charge_cutoff + 1);
nothing #hide
```

These circuit parameters set the capacitive coupling coefficients. `op` names
a component operator and `param` names a model parameter; the component names
supply the prefixes `q1`, `q2`, and `c`.

```@example coupler
η = EC12 * ECc / (EC1c * EC2c)
g1c = 8 * EC1 * ECc / EC1c
g2c = 8 * EC2 * ECc / EC2c
g12 = 8 * (1 + η) * EC1 * EC2 / EC12

interaction = param(:g1c) * op(:q1_charge) * op(:c_charge) +
    param(:g2c) * op(:q2_charge) * op(:c_charge) +
    param(:g12) * op(:q1_charge) * op(:q2_charge)
model = make_model([q1, q2, coupler], interaction, (; g1c, g2c, g12);
    truncation_dimensions=Dict(q1 => retained_levels, q2 => retained_levels,
        coupler => retained_levels),
    max_dimension=retained_levels^3)
size(model.H)
```

Each transmon starts in a 121-state charge basis and retains three levels,
giving a 27-dimensional coupled model. This keeps the example small; increase
the retained levels and check convergence before using the pulse quantitatively.
Pulsed flux changes the Hamiltonian in the basis constructed at idle flux.

`model.states` contains dressed eigenstates of the coupled idle device. Labels
are zero-based and follow the component order `[q1, q2, coupler]`.

| Label | Dressed state |
|:------|:--------------|
| `(0,0,0)` | Both qubits and coupler in the ground state |
| `(1,0,0)` | One excitation in qubit 1 |
| `(0,1,0)` | One excitation in qubit 2 |
| `(1,1,0)` | Both qubits excited, coupler in the ground state |

```@example coupler
exchange_states = [model.states[(1, 0, 0)], model.states[(0, 1, 0)]];
nothing #hide
```

## 3. Define the flux pulse

Qubit 1 is higher in frequency than qubit 2 at idle. Its flux excursion lowers
its frequency, while the coupler excursion changes the mediated interaction.
These amplitudes are starting guesses near a working exchange pulse.
We keep the duration and sin² ramp times fixed during calibration.

```@example coupler
duration = 28.135660302912168
ramp = 1.8727503955158111
q1_flux = Pulse(RampedFlattop(ramp); amplitude=0.130, duration)
coupler_flux = Pulse(RampedFlattop(ramp); amplitude=0.385,
    delay=ramp, duration=duration - 2ramp)
model.gates[:exchange] = DeviceGate((; q1_phi=q1_flux, c_phi=coupler_flux), 0);
nothing #hide
```

`q1_phi` and `c_phi` replace existing model parameters during evolution, so the
gate needs zero additional Hamiltonian. Both pulses return to the idle value
of zero at the endpoints. The coupler starts one ramp later and finishes one
ramp earlier than qubit 1.

Plot the controls to see this timing before simulating the gate.

```@example coupler
times = range(0, duration; length=201)
fluxes = Dict(:q1 => q1_flux.(times), :coupler => coupler_flux.(times))
flux_styles = Dict(:q1 => (; color=:royalblue, label="Qubit 1"),
    :coupler => (; color=:purple, label="Coupler"))
flux_figure, ax = plot_trajectories(fluxes, times, flux_styles;
    axis=(; xlabel="Time (ns)", ylabel="Flux (Φ₀)", title="Trial flux pulse"),
    legend_options=(; position=:lt))
flux_figure
```

## 4. Simulate and project the exchange

Start with one excitation in qubit 1. Evolve the full model and measure the
amplitude in the dressed target state; its squared magnitude is the population.

```@example coupler
evolution_options = (; abstol=1e-8, reltol=1e-8)
result = sesolve(model, model.gates[:exchange], exchange_states[1], times;
    progress_bar=false, evolution_options...)
amplitudes = state_amplitudes(Dict(:target => exchange_states[2]), result)
transfer_before = abs2.(amplitudes[:target])
last(transfer_before)
```

To check both transfer directions, supply the ordered pair of states to
`numerical`. It evolves both input states through the full retained model and
then projects onto this pair. The result is a 2×2 operator; `.data` gives its
matrix. Projection does not renormalize columns, so leakage remains visible.

```@example coupler
U_exchange = numerical(model, model.gates[:exchange], exchange_states;
    frame=:lab, evolution_options...)
U_exchange.data
```

Within this pair, an X matrix exchanges the two populations. With
`include_phases=false`, the score is
$(|U_{21}|^2+|U_{12}|^2)/2$: the mean transfer probability over both inputs.
It ignores relative phases and is therefore a population score rather than a
coherent gate fidelity. Keep the target `sigmax()` unscaled.

```@example coupler
unitary_fidelity(sigmax(), U_exchange; include_phases=false)
```

## 5. Calibrate the flux amplitudes

Select the two amplitudes by their parameter paths. `CalibrationProblem`
prepares the loss, `solve` runs Nelder–Mead, and `calibrated_gate` reconstructs
the pulse. A small initial simplex suits these nearby flux guesses.
The objective uses both exchange directions and includes leakage penalties.

```@example coupler
selectors = ["q1_phi/amplitude", "c_phi/amplitude"]
simplex = OptimizationOptimJL.Optim.AffineSimplexer(a=0.001, b=0.001)
algorithm = NelderMead(initial_simplex=simplex)
problem = CalibrationProblem(model, model.gates[:exchange], selectors, sigmax();
    states=exchange_states, frame=:lab, include_phases=false,
    evolution_kwargs=evolution_options)
solution = solve(problem, algorithm; maxiters=100, g_abstol=1e-8)
model.gates[:exchange] = calibrated_gate(problem, solution.u)

U_exchange = numerical(model, model.gates[:exchange], exchange_states;
    frame=:lab, evolution_options...)
unitary_fidelity(sigmax(), U_exchange; include_phases=false)
```

Inspect the returned score and `solution.retcode`, then compare the target
population before and after calibration. With these settings, calibration
raises the transfer probability from approximately 0.844 to 0.996.

```@example coupler
result = sesolve(model, model.gates[:exchange], exchange_states[1], times;
    progress_bar=false, evolution_options...)
amplitudes = state_amplitudes(Dict(:target => exchange_states[2]), result)
transfer_after = abs2.(amplitudes[:target])

styles = Dict(:trial => (; color=:gray, label="Trial"),
    :calibrated => (; color=:purple, label="Calibrated"))
transfer_figure, ax = plot_trajectories(
    Dict(:trial => transfer_before, :calibrated => transfer_after), times, styles;
    axis=(; xlabel="Time (ns)", ylabel="Target population", title="Coupler exchange"),
    legend_options=(; position=:lt))
transfer_figure
```

### Check the full computational subspace

Good excitation transfer alone does not establish a two-qubit SWAP. Include
the states with zero and two qubit excitations, keeping the coupler in its
ground state, and compare the resulting 4×4 matrix with a full SWAP target.
The basis order below is $|00\rangle, |01\rangle, |10\rangle, |11\rangle$.

```@example coupler
computational_states = [model.states[(0, 0, 0)], model.states[(0, 1, 0)],
    model.states[(1, 0, 0)], model.states[(1, 1, 0)]]
U_computational = numerical(model, model.gates[:exchange], computational_states;
    frame=:lab, evolution_options...)
U_computational.data

swap_target = QuantumObject([1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1])
unitary_fidelity(swap_target, U_computational)
```

This coherent score retains relative phases and ignores only a common global
phase. Our calibration optimized population exchange, so it need not give a
high SWAP fidelity: the score here is approximately 0.323. A complete gate
calibration must also address phases and
the action on $|11\rangle$. Use `frame=:interaction` consistently in calibration
and evaluation if you want to remove idle evolution; that alone does not
remove phases acquired during the flux excursion.

## 6. Add a microwave control

Microwave pulses add a charge-drive term instead of replacing a flux parameter.
Here is a starting pulse for a π/2 rotation on qubit 1. Its carrier frequency
comes from the dressed idle spectrum, and `SineSquared` gives a single smooth
lobe over the pulse duration.

```@example coupler
q1_frequency = abs(model.spectrum[(1, 0, 0)] - model.spectrum[(0, 0, 0)])
drive = Pulse(SineSquared(); amplitude=0.024775449345990953,
    duration=9.027992059968785,
    carrier=SineCarrier(q1_frequency; phase=1.5501484428191823))
model.gates[:quarter_x1] = DeviceGate((; drive), param(:drive) * op(:q1_charge));
nothing #hide
```

To calibrate this pulse, select its amplitude, frequency, and carrier phase as
in [Transmon resonator control](transmon_resonator_control.md). For a two-qubit
single-qubit gate, check both spectator-qubit states in the computational
subspace. The simulations here are closed-system and omit dissipation.

The model now contains the calibrated flux gate and the microwave pulse. See
[Saving and loading model bundles](../user_guide/persistence.md) to reuse them.
