"""
    AbstractCalibrationSetup

Abstract interface for calibration setup objects. The built-in implementation
[`CalibrationProblem`](@ref) exposes a standard SciML optimization problem and
a callable that reconstructs gates from candidate values.
"""
abstract type AbstractCalibrationSetup end

"""
    CalibrationProblem(problem, rebuild)
    CalibrationProblem(model, gate, names, target; objective=gate_infidelity,
        states=nothing, output_states=states, frame=:lab, include_phases=true, energy_shift=:auto, dense=false,
        evolution_kwargs=(;), kwargs...)

Calibration problem containing a standard `SciMLBase.OptimizationProblem`
and `rebuild(values) -> DeviceGate`. Select parameters with
`CalibrationProblem(model, gate, names, target; kwargs...)`, then call
`solve(problem, algorithm; kwargs...)`. This forwards directly to SciML and returns its optimization solution unchanged. Reconstruct the gate with
[`calibrated_gate`](@ref).
Algorithm choice and solver options remain the caller's responsibility.
"""
struct CalibrationProblem{P,F} <: AbstractCalibrationSetup
    problem::P
    rebuild::F
end

"""
    SciMLCalibrationSetup

Compatibility alias for [`CalibrationProblem`](@ref).
"""
const SciMLCalibrationSetup = CalibrationProblem

"""
    parameters(object)

List named parameters as `name => (path, value)`. Gate scalar names are unchanged;
pulse parameters use slash-separated names such as `"drive/amplitude"`.
Values such as `nothing` are shown, but calibration only accepts finite real values.
"""
function parameters end

function _parameter_name(name)
    name = String(name)
    isempty(name) || occursin('/', name) ?
        throw(ArgumentError("Parameter names must be nonempty and cannot contain '/'")) : name
end

# Explicitly traverse .parameters when a parameter shadows a structure field.
_parameter_path(object, key) = hasfield(typeof(object), key) ?
    "parameters/$(key)" : String(key)

# Kept as a discovery alias; public calibration paths are now strings.
"""
    calibration_values(object)

Discovery alias for [`parameters`](@ref). Returns the same mapping of parameter
names to `(path, value)` tuples; it does not filter values to optimizable scalars.
Calibration accepts only selected finite real values.
"""
calibration_values(object) = parameters(object)

function _selected_calibration_parameters(gate::DeviceGate, names)
    available = parameters(gate)
    selected_values = Float64[]
    paths = String[]
    for name in names
        name isa AbstractString || throw(ArgumentError("Calibration parameter names must be strings"))
        haskey(available, name) || throw(ArgumentError("Unknown calibration parameter: $name"))
        path, value = available[name]
        value isa Real && isfinite(value) || throw(ArgumentError("Calibration parameter $name must be a finite real scalar"))
        isequal(getpath(gate, path), value) || throw(ArgumentError("Parameter $name does not match its path $path"))
        push!(selected_values, value)
        push!(paths, path)
    end
    return selected_values, paths
end

function _with_calibration_values(gate::DeviceGate, paths, values)
    length(values) == length(paths) || throw(DimensionMismatch(
        "Expected $(length(paths)) calibration values, received $(length(values))",
    ))
    return foldl(zip(paths, values); init = gate) do candidate, (path, value)
        setpath(candidate, path, value)
    end
end

"""
    CalibrationProblem(model, gate, names, target; objective=gate_infidelity,
        states=nothing, output_states=states, frame=:lab, include_phases=true, energy_shift=:auto, dense=false,
        evolution_kwargs=(;), kwargs...)

Build a [`CalibrationProblem`](@ref) wrapping a standard
`SciMLBase.OptimizationProblem` and an immutable gate rebuilder.
Select parameters using names from `parameters(gate)`, e.g. `["drive/amplitude"]`.
Additional keyword arguments (including `lb` and `ub`) are forwarded to
`OptimizationProblem`. The objective has signature `(model, candidate_gate, target)`
and must return a scalar loss. The default compares full retained-space unitaries.
Supply `states` to compare a selected gate subspace; `output_states` defaults to
that basis. `frame=:interaction` removes idle evolution before comparison.
`include_phases=false` compares probabilities instead of coherent gate matrices.
`energy_shift=:auto` caches the mean idle energy of the input states, subtracts
it during integration, and restores global phase before gate comparison.
An explicit offset uses frequency units; `energy_shift=0` disables centering.
`dense=true` constructs dense Hamiltonian matrices, including the cached fixed terms.
The default `dense=false` preserves their usual storage.
`evolution_kwargs` forwards native evolution options. These built-in options
cannot be combined with a custom objective. Operators, fixed terms, bases, and
target validation are prepared once per problem, preserving sparse storage.

Selected names must be a nonempty collection of unique discovery names with finite
real values. The input gate is unchanged. Solve with
`solve(problem, algorithm; kwargs...)`; solver options are forwarded to SciML.
"""
function CalibrationProblem(model::DeviceModel, gate::DeviceGate, names, target;
    objective=gate_infidelity, states=nothing, output_states=states, frame=:lab,
    include_phases::Bool=true, energy_shift=:auto, dense::Bool=false, evolution_kwargs=(;), kwargs...)
    names isa AbstractString && throw(ArgumentError("Pass a list of parameter names, such as [\"drive/amplitude\"]"))
    names = collect(names)
    isempty(names) && throw(ArgumentError("At least one calibration parameter is required"))
    allunique(names) || throw(ArgumentError("Calibration parameter names must be unique"))
    initial_values, paths = _selected_calibration_parameters(gate, names)
    rebuild = values -> _with_calibration_values(gate, paths, values)
    if objective === gate_infidelity
        prepared = _prepared_gate_objective(model, gate, target;
            states, output_states, frame, include_phases, energy_shift, dense, evolution_kwargs)
        loss = (values, _) -> prepared(rebuild(values))
    else
        (states === nothing && output_states === nothing && frame === :lab &&
            include_phases && energy_shift === :auto && !dense && isempty(evolution_kwargs)) || throw(ArgumentError(
            "Custom objectives cannot be combined with built-in gate-evolution options"))
        loss = (values, _) -> objective(model, rebuild(values), target)
    end
    problem = SciMLBase.OptimizationProblem(loss, initial_values; kwargs...)
    return CalibrationProblem(problem, rebuild)
end

"""
    calibration_problem(model, gate, names, target; kwargs...)

Compatibility constructor for [`CalibrationProblem`](@ref). All arguments and
keywords are passed to that constructor unchanged.
"""
calibration_problem(model::DeviceModel, gate::DeviceGate, names, target; kwargs...) =
    CalibrationProblem(model, gate, names, target; kwargs...)

"""
    calibrated_gate(problem::CalibrationProblem, values)

Reconstruct a gate from candidate parameter values in the selected order without
modifying the original gate. The vector length must match the calibration problem.
"""
calibrated_gate(problem::CalibrationProblem, values) = problem.rebuild(values)

"""
    calibrate(setup, algorithm; kwargs...)

Solve a calibration through the standard SciML interface and rebuild its gate.
All keyword arguments are forwarded directly to `SciMLBase.solve`.
Returns `(solution, gate)` with the candidate reconstructed from `solution.u`.
Inspect the solution retcode and objective before accepting the result; this wrapper
does not independently enforce optimization convergence. The original gate is
unchanged. Prefer `solution = solve(problem, algorithm; kwargs...)` followed by
`calibrated_gate(problem, solution.u)`. This function is retained for compatibility.
"""
function calibrate(setup::CalibrationProblem, algorithm; kwargs...)
    solution = SciMLBase.solve(setup, algorithm; kwargs...)
    return solution, calibrated_gate(setup, solution.u)
end
