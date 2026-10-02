# Tunable coupler control

Build two data transmons and a flux-tunable coupler, define flux and microwave
controls, and simulate excitation exchange between the qubits. This tutorial
uses fixed circuit and pulse parameters to demonstrate the package workflow.

## 1. Set up the tutorial

From the repository root, install the demo dependencies once:

```sh
julia --project=demo -e 'using Pkg; Pkg.instantiate()'
```

Run the complete [Julia tutorial](https://github.com/Gavin-Rockwood/QuantumDevices.jl/blob/main/demo/tunable_coupler.jl)
with `julia --project=demo demo/tunable_coupler.jl`. It prints the final target
population and saves a model bundle and two figures in a temporary directory.
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
cycles/ns (GHz), while time is in ns. Keep these values in frequency units
when constructing Hamiltonians; apply `2π` only at evolution.

## 2. Construct and couple three transmons

The model helper constructs `q1`, `q2`, and `c` with
`make_tunable_transmon`, then couples their charge operators using

```julia
interaction = param(:g1c) * op(:q1_charge) * op(:c_charge) +
              param(:g2c) * op(:q2_charge) * op(:c_charge) +
              param(:g12) * op(:q1_charge) * op(:q2_charge)
```

Operator prefixes come from the component names. The helper computes coupling
strengths from `TARGET` in GHz and passes the interaction
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

`make_tunable_transmon` takes the charging energy and the two junction
energies in frequency units (`E/h`, GHz here). This model splits each total Josephson energy equally
between the junctions. Its `phi` parameter is flux in units of the flux quantum;
use that same parameter directly when defining controls.

## 3. Drive model parameters with flux pulses

A flux control replaces an existing model parameter over the gate duration.
Build the exchange pulse explicitly:

```@example coupler
settings = SWAP_SETTINGS
T, ramp = settings.duration, settings.ramp
q1_flux = Pulse(RampedFlattop(ramp); amplitude=settings.q1_flux, duration=T)
c_flux = Pulse(RampedFlattop(ramp); amplitude=settings.coupler_flux,
    delay=ramp, duration=T-2ramp)
flux_gate = DeviceGate((; q1_phi=q1_flux, c_phi=c_flux), 0)
@assert flux_gate.duration == T

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
   midpoint_coefficient=x1.parameters.drive(x1.duration / 2))
```

Inside `quarter_gate`, `Pulse(SineSquared(); duration=T, amplitude=epsilon,
carrier=SineCarrier(frequency; phase))` defines a microwave drive directly.
The envelope has one lobe over its duration; carrier phase is in radians.


## 5. Simulate excitation exchange

Start in the dressed `|100⟩` state and evolve under the flux gate:

```@example coupler
gate = model.gates[:swap]
source = model.states[(1, 0, 0)]
target = model.states[(0, 1, 0)]
times = range(0, gate.duration; length=121)
result = sesolve(2pi * numerical(model, gate), source, times;
    progress_bar=false, tstops=pulse_tstops(gate), abstol=1e-8, reltol=1e-8)
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
simulation. The `:swap` gate transfers an excitation. Establishing a full SWAP
gate also requires checking phases and the action on the other computational
states. This tutorial measures transfer without adding dissipation or estimating
a gate fidelity.

## 6. Save and reload the model

Save the instantiated model with its flux gate and four microwave gates:

```@example coupler
bundle_path = QuantumDevices.save(joinpath(mktempdir(), "tunable_coupler"), model)
restored = QuantumDevices.load(bundle_path)
@assert restored.H ≈ model.H
@assert Set(keys(restored.gates)) == Set(keys(model.gates))
@assert restored.gates[:quarter_x1].parameters.drive(x1.duration / 2) ≈
    x1.parameters.drive(x1.duration / 2)
println(join(sort(readdir(bundle_path)), ", "))
```

`QuantumDevices.save` creates a new directory and refuses an existing destination.
The standalone demo writes `<output_directory>/model/` alongside `controls.png`
and `transfer.png`, and prints the bundle path. All circuit and pulse parameters
are defined directly in `demo/tunable_coupler.jl`.

To reload in a fresh Julia session, include the module defining its charge-drive
callable first:

```julia
using QuantumDevices
include(joinpath(pkgdir(QuantumDevices), "demo", "tunable_coupler.jl"))
restored = QuantumDevices.load("output_directory/model")
```

Including the file defines its module without running or saving the tutorial.
See [Saving and loading model bundles](../user_guide/persistence.md) for the
bundle layout.
