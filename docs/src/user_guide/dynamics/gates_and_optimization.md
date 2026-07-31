# Gates and Optimization

`GateSpec` binds finite-duration trajectories to the non-fixed controls exposed
by one `ModelSpec`. A gate starts and ends at the model defaults, so it can be
composed with idle evolution without a Hamiltonian discontinuity.

```julia
drive_path = ParamPath(:drive)
gate = GateSpec(
    :x_gate,
    spec;
    duration = 10.0,
    parameters = (amplitude = 0.1, sigma = 2.0),
    controls = Dict(
        drive_path => CarrierControl(
            (p, t) -> p.amplitude * gaussian(t, p.duration, p.sigma);
            frequency = 4.8,
        ),
    ),
)

H_gate = numerical(gate)
stronger_gate = with_parameters(gate; amplitude = 0.12)
```

Only paths returned by `parameters(model(spec))` may be controlled. Fixed,
unknown, out-of-domain, and non-idle endpoint controls are rejected.

Registering a gate first verifies that its `ModelSpec` is the exact recipe
registered on the device. The stored gate then receives an independent snapshot
of that recipe, so later component or model updates cannot change it.

Parameterized controls receive `p.duration` and the gate's named scalar
parameters. Existing `t -> value` controls remain supported. `with_parameters`
reconstructs a gate without mutating the original.

`CarrierControl` keeps a slow analytic envelope separate from its fixed
laboratory-frame carrier. With frequencies in GHz and time in ns, its resolved
coefficient is `envelope * cos(2π * frequency * t + phase)`. Classical
propagation sees the complete waveform, while the Piccolo extension optimizes
the slow envelope through a time-modulated drive.

## Gate optimization

The common workflow optimizes named recipe parameters directly:

```julia
result = optimize(
    gate,
    ComplexF64[0 1; 1 0],
    :amplitude => (0.0, 0.2);
    states = [(0,), (1,)],
    frame_frequencies = (q = 4.8,),
    iterations = 80,
)
```

Multiple selectors use a tuple of pairs. Bounds may instead be passed
separately as a dictionary or named tuple. `method` accepts Optim algorithms,
`options` accepts `Optim.Options`, and propagation settings belong in
`solver_kwargs`.

The default backend evaluates deterministic coarse seeds and refines the best
one with bounded derivative-free optimization. The result retains the optimized
gate recipe, named minimizer, fidelity, leakage, failures, settings, and backend
result. The input gate is never mutated.

`StateTransferObjective(initial, target)` uses final-state overlap. States and
logical subspaces may be specified with the ordered state labels owned by the
model. For lab-frame unitary propagation, `frame_frequencies` applies the known
idle-frame virtual-Z correction before process fidelity is scored. Its named
keys must match the model subsystem order.

### Low-level problem API

An `EnvelopeOptimizationProblem` uses named bounded scalar variables and a
builder that produces candidate gates. It remains available for custom builders:

```julia
variables = [
    GateVariable(:amplitude; initial = 0.2, lower = 0.0, upper = 1.0),
]

objective = UnitaryObjective(
    ComplexF64[0 1; 1 0];
    states = [(0,), (1,)],
)

problem = EnvelopeOptimizationProblem(
    values -> GateSpec(
        :x_gate,
        spec;
        duration = 1.0,
        controls = Dict(
            drive_path =>
                (t -> values.amplitude * sinpi(t)^2),
        ),
    ),
    variables,
    objective,
)

result = optimize_gate(problem)
```

## Piccolo

Piccolo support loads only when `Piccolo.jl` is installed and loaded. Select
model controls by `ParamPath`:

```julia
using Piccolo

result = optimize(
    gate,
    ComplexF64[0 1; 1 0],
    drive_path => (-1.0, 1.0);
    backend = :piccolo,
    states = [(0,), (1,)],
    knots = 50,
)
```

The extension converts the symbolic parameterized model into Hermitian drift
and linear or nonlinear Piccolo drive terms, optimizes a smooth
piecewise-constant pulse, and converts the result into a resolved `GateSpec`.
Add `:duration => (lower, upper)` to run Piccolo's minimum-time wrapper after
the smooth-pulse solve. `PiccoloOptimizationProblem` and `optimize_gate` remain
available as the low-level interface.

Carrier-modulated Piccolo controls require zero model defaults and linear
control terms. Minimum-time optimization is rejected for frame-corrected
unitary targets because their corresponding lab-frame target changes with gate
duration.

See the
[physical gate optimization notebook](https://github.com/Gavin-Rockwood/QuantumDevices.jl/blob/main/demo/GateOptimizationExamples.ipynb)
for a complete device-first workflow: one physical lab-frame model, Gaussian
X-gate and cos²-ramped iSWAP calibration, gate registration, independent
QuantumToolbox verification, and optional carrier-aware Piccolo synthesis.
