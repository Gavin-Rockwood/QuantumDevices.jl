# Tunable coupler control

Build two data transmons and a flux-tunable coupler, define flux and microwave
controls, and simulate excitation exchange between the qubits. This tutorial
uses fixed circuit and pulse parameters from an existing simulation.

## 1. Set up the demo

From the repository root, install the demo dependencies once:

```sh
julia --project=demo -e 'using Pkg; Pkg.instantiate()'
```

Run the complete [Julia tutorial](https://github.com/Gavin-Rockwood/QuantumDevices.jl/blob/main/demo/tunable_coupler.jl)
with `julia --project=demo demo/tunable_coupler.jl`. It prints the final target
population and saves control and transfer figures in a temporary directory.
Pass a directory as the last argument to choose the output location. For an
interactive session, start `julia --project=demo` and follow the steps below.

```@example coupler
using QuantumDevices, QuantumToolbox, CairoMakie, LinearAlgebra, SciMLBase
include(joinpath(pkgdir(QuantumDevices), "demo", "tunable_coupler.jl"))
using .TunableCouplerDemo
TARGET
```

`TARGET` contains circuit parameters; `SWAP_SETTINGS` and `QUARTER_SETTINGS`
contain pulse parameters. Energies, couplings, and drive frequencies are in
cycles/ns (GHz), while time is in ns. Convert energy and drive coefficients to
angular units with `2π` when constructing Hamiltonians.

## 2. Construct and couple three transmons

The model helper constructs `q1`, `q2`, and `c` with
`make_tunable_transmon`, then couples their charge operators using

```julia
interaction = param(:g1c) * op(:q1_charge) * op(:c_charge) +
              param(:g2c) * op(:q2_charge) * op(:c_charge) +
              param(:g12) * op(:q1_charge) * op(:q2_charge)
```

Operator prefixes come from the component names. The helper computes coupling
strengths from `TARGET`, converts them to rad/ns, and passes the interaction
and parameter values to `make_model`.

```@example coupler
model = make_coupler_model(; retained_levels=3, n_cutoff=60)
@assert size(model.H) == (27, 27)
(; hamiltonian_size=size(model.H), gates=collect(keys(model.gates)))
```

Each transmon starts in a 121-state charge basis and retains three levels,
giving a 27-dimensional coupled space. Operators such as `n²` are projected
from the parent basis before truncation. The helper also attaches one flux
pulse and four microwave pulses.

The original flux convention uses
`EJ_eff(Φ) = EJ sqrt(cos(2πΦ)^2 + d^2 sin(2πΦ)^2)`.
`make_tunable_transmon` uses `cos(π phi)`, so these controls set `phi = 2Φ`.
Equal junction energies give the original `d=0`: each constructor receives
`EJ1 = EJ2 = π EJ` in angular units.

## 3. Drive model parameters with flux pulses

A flux control replaces an existing model parameter over the gate duration.
Build the exchange pulse explicitly:

```@example coupler
settings = SWAP_SETTINGS
T, ramp = settings.duration, settings.ramp
q1_flux = ramped_flattop_pulse(2 * settings.q1_flux, ramp;
    ramp=:linear, stop=T)
c_flux = ramped_flattop_pulse(2 * settings.coupler_flux, ramp;
    ramp=:linear, start=ramp, stop=T-ramp)
model.gates[:swap] = DeviceGate((; q1_phi=q1_flux, c_phi=c_flux), 0, T)
model.gates[:swap].duration
```

`q1_phi` and `c_phi` match parameters already present in the model Hamiltonian.
The gate therefore needs zero additional Hamiltonian. The coupler starts one
ramp later than qubit 1 and finishes one ramp earlier. Both pulses return to
the idle flux at the gate endpoints. The `swap_gate()` helper constructs the
same gate.

## 4. Add microwave controls

Microwave gates add a time-dependent charge drive. For example, select the
quarter-X drive on qubit 1:

```@example coupler
x1 = quarter_gate(1, :X)
model.gates[:quarter_x1] = x1
(; duration_ns=x1.duration,
   midpoint_coefficient=pulse_value(x1.parameters.drive, x1.duration / 2, x1.duration))
```

Inside `quarter_gate`, `sine_squared_pulse(1.0, 1/(2T))` gives one envelope lobe
over the gate duration. A `GenericPulseFunction` combines it with a carrier,

```math
h(t)=2\pi\epsilon\,e(t)\sin(2\pi\nu t+\varphi),
```

and `DeviceGate((; drive), param(:drive) * op(:q1_charge), T)` attaches that
coefficient to qubit 1. The callable receives the owning gate duration through
`p.duration`. Choose qubit `2` or axis `:Y` to build the other controls. Their
names describe the intended rotations from the original parameters; this
tutorial does not recalibrate or verify those rotations.

```@example coupler
plot_controls(model)
```

The upper panel shows flux in the original `Φ` convention, dividing model
`phi` by two. The lower panels show microwave coefficients divided by `2π`.

## 5. Simulate excitation exchange

Start in the dressed `|100⟩` state and evolve under the flux gate:

```@example coupler
gate = model.gates[:swap]
source = model.states[(1, 0, 0)]
target = model.states[(0, 1, 0)]
times = range(0, gate.duration; length=121)
result = sesolve(numerical(model, gate), source, times;
    progress_bar=false, abstol=1e-8, reltol=1e-8)
@assert SciMLBase.successful_retcode(result.retcode)
p100 = [abs2(dot(source, state)) for state in result.states]
p010 = [abs2(dot(target, state)) for state in result.states]
@assert last(p010) > 0.99
(; final_source_population=last(p100), final_target_population=last(p010))
```

State labels are zero-based and follow `[q1, q2, c]`. `numerical(model, gate)`
combines the model Hamiltonian with its time-dependent controls for `sesolve`.
Squared overlaps measure the source and target populations. The script's
`swap_population(model)` helper performs these same steps.

```@example coupler
transfer = Figure(size=(700, 370))
ax = Axis(transfer[1, 1]; xlabel="Time (ns)", ylabel="Dressed-state population",
    title="Coupler flux pulse")
lines!(ax, times, p100; label="|100⟩", color=:royalblue, linewidth=2.5)
lines!(ax, times, p010; label="|010⟩", color=:purple, linewidth=2.5)
axislegend(ax; position=:rt)
ylims!(ax, 0, 1.05)
transfer
```

You should obtain a final `|010⟩` population above 0.99 in this closed-system
simulation. The original pulse is labeled `SWAP`, but population transfer alone
does not establish a full SWAP gate: phases and the action on the other
computational states also matter. This tutorial measures transfer without
adding dissipation or estimating a gate fidelity.
