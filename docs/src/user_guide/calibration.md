# Calibrating a gate

Calibration selects finite real parameters, evaluates an objective, and rebuilds
a new gate. QuantumDevices supplies the standard SciML optimization problem;
you choose the optimizer and its options. Install OptimizationOptimJL in your
application environment for this example.

```@example calibration
using QuantumDevices, QuantumToolbox, SciMLBase, OptimizationOptimJL, CairoMakie
q = make_qubit("q", 0.0)
model = make_model([q], val(0), (;))
trial = DeviceGate((; drive=constant_pulse(0.4)),
                   param(:drive) * op(:q_x), 1.0)
target = -im * sigmax()
@assert haskey(parameters(trial), "drive/amplitude")
parameters(trial)["drive/amplitude"]
```

The ideal zero-idle Hamiltonian is `A X`, so `U(T)=exp(-im*A*T*X)`.
The target is reached at `A*T=π/2` up to global phase.

```@example calibration
errors = Float64[]
objective = function (m, gate, target)
    value = gate_infidelity(m, gate, target)
    push!(errors, value)
    value
end
setup = calibration_problem(model, trial, ["drive/amplitude"], target; objective)
solution, calibrated = calibrate(setup, OptimizationOptimJL.NelderMead(); maxiters=80)
@assert SciMLBase.successful_retcode(solution.retcode)
@assert solution.objective < 1e-7
@assert isapprox(getpath(calibrated, "drive/amplitude"), π/2; atol=1e-3)
@assert getpath(trial, "drive/amplitude") == 0.4
getpath(calibrated, "drive/amplitude")
```

Objective kwargs such as `lb` and `ub` go to `calibration_problem`; optimizer
options such as `maxiters` go to `calibrate`. Select an algorithm that supports
your bounds and differentiation strategy. Non-real parameters and undiscovered
names cannot be optimized by this wrapper.

```@example calibration
fig = Figure(size=(700, 320))
ax = Axis(fig[1, 1], xlabel="Objective evaluation", ylabel="Best process infidelity",
          yscale=log10)
lines!(ax, eachindex(errors), max.(accumulate(min, errors), 1e-14))
fig
```

The plotted floor is for display only; assertions use the computed objective.

## Direct SciML access and custom objectives

`setup.problem` can be passed directly to `SciMLBase.solve`. Reconstruct a gate
with `calibrated_gate(setup, solution.u)`. This produces the same result as the
wrapper without a second optimization.

Custom objectives have signature `objective(model, candidate_gate, target)` and
return a scalar loss. They may represent a subspace target or leakage penalty,
but their definition and validation are your responsibility. The built-in
`gate_infidelity` compares full retained-space unitaries, so it includes all levels
kept by the model. See [fidelity conventions](../explanations/conventions.md).
