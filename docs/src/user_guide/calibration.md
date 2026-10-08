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
states = [basis(2, 0), basis(2, 1)]
problem = CalibrationProblem(model, trial, ["drive/amplitude"], target; states, dense=true)
solution = solve(problem, OptimizationOptimJL.NelderMead(); maxiters=80, store_trace=true)
errors = OptimizationOptimJL.Optim.f_trace(solution.original)
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
ax = Axis(fig[1, 1], xlabel="Optimizer iteration", ylabel="Best process infidelity",
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
but their definition and validation are your responsibility. The built-in objective evolves the full retained space by default.
Supply `states=[ket0, ket1, ...]` to evolve just those ordered columns and compare
a target matrix in that basis. `output_states` defaults to the same basis.

`include_phases=true` includes relative phases and ignores global phase. Setting
it to `false` matches transition probabilities; leakage still lowers the score.
The default `frame=:lab` preserves lab phases. Explicit `frame=:interaction`
removes idle evolution from the final comparison and leaves integration unchanged.
Pass native evolution tolerances through `evolution_kwargs=(; abstol=1e-8, reltol=1e-8)`. See [fidelity conventions](../explanations/conventions.md).

Calibration automatically subtracts the mean idle energy of the selected input
states during integration. The offset is cached when `CalibrationProblem` is
constructed, and global phase is restored before comparison. This changes neither
the lab-frame target nor the `frame=:interaction` convention. Use `energy_shift=0`
to disable centering, or supply a fixed offset in frequency units.

Set `dense=true` directly on `CalibrationProblem` to use dense Hamiltonian matrices
for every candidate. This includes the fixed matrices cached when the problem is
constructed. The default `dense=false` preserves the usual storage. Dense storage
can be faster for small Hamiltonians with relatively few zeros; benchmark the
choice on your model. The `solve` call and native SciML options stay the same.
