# Calibrating a gate

Calibration selects finite real parameters, evaluates an objective, and rebuilds
a new gate. QuantumDevices supplies the standard SciML optimization problem;
you choose the optimizer and its options. Install OptimizationOptimJL in your
application environment for this example.

```@example calibration
using QuantumDevices, QuantumToolbox, SciMLBase, OptimizationOptimJL, CairoMakie
q = make_qubit("q", 0.0)
model = make_model([q], val(0), (;))
trial = DeviceGate((; drive=Pulse(Constant(); amplitude=0.1, duration=1.0)),
                   param(:drive) * op(:q_x), 1.0)
target = -im * sigmax()
@assert haskey(parameters(trial), "drive/amplitude")
parameters(trial)["drive/amplitude"]
```

The ideal zero-idle Hamiltonian is `A X`, so `U(T)=exp(-2π*im*A*T*X)`.
The target is reached at `A*T=1/4` up to global phase.

```@example calibration
errors = Float64[]
objective = function (m, gate, target)
    value = gate_infidelity(m, gate, target)
    push!(errors, value)
    value
end
problem = CalibrationProblem(model, trial, ["drive/amplitude"], target; objective)
solution = solve(problem, OptimizationOptimJL.NelderMead(); maxiters=80)
calibrated = calibrated_gate(problem, solution.u)
@assert SciMLBase.successful_retcode(solution.retcode)
@assert solution.objective < 1e-7
@assert isapprox(getpath(calibrated, "drive/amplitude"), 0.25; atol=1e-3)
@assert getpath(trial, "drive/amplitude") == 0.1
getpath(calibrated, "drive/amplitude")
```

Problem options such as `lb` and `ub` go to `CalibrationProblem`; optimizer
options such as `maxiters` go to `solve`. Select an algorithm that supports
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

`solve(problem, algorithm; kwargs...)` forwards directly to SciML and returns its
standard optimization solution. All solver options, including callbacks and
tolerances, pass through unchanged. `problem.problem` exposes the underlying
`OptimizationProblem` when needed. `calibration_problem(...)` and `calibrate(...)`
remain available for compatibility.

Custom objectives have signature `objective(model, candidate_gate, target)` and
return a scalar loss. They may represent a subspace target or leakage penalty,
but their definition and validation are your responsibility. The built-in
`gate_infidelity` compares full retained-space unitaries, so it includes all levels
kept by the model. See [fidelity conventions](../explanations/conventions.md).
