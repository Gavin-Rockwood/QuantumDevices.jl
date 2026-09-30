"""
    AbstractCalibrationSetup

Abstract interface for calibration setup objects. The built-in implementation
[`SciMLCalibrationSetup`](@ref) exposes a standard SciML optimization problem and
a callable that reconstructs gates from candidate values.
"""
abstract type AbstractCalibrationSetup end

"""
    SciMLCalibrationSetup(problem, rebuild)

Calibration result setup containing a standard `SciMLBase.OptimizationProblem`
and `rebuild(values) -> DeviceGate`. Construct it with [`calibration_problem`](@ref),
then call [`calibrate`](@ref) or `SciMLBase.solve(setup.problem, algorithm)`.
Algorithm choice and solver options remain the caller's responsibility.
"""
struct SciMLCalibrationSetup{P,F} <: AbstractCalibrationSetup
    problem::P
    rebuild::F
end

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

parameters(pulse::AbstractParameterizedPulse) = Dict{String,Tuple}(
    _parameter_name(name) => (_parameter_path(pulse, name), value)
    for (name, value) in pairs(pulse.parameters))
parameters(::AbstractPulse) = Dict{String,Tuple}()

function parameters(gate::DeviceGate)
    result = Dict{String,Tuple}()
    for (name, value) in pairs(gate.parameters)
        name_string = _parameter_name(name)
        path = _parameter_path(gate, name)
        if value isa AbstractPulse
            for (pulse_name, (pulse_path, pulse_value)) in parameters(value)
                full_path = "$path/$pulse_path"
                haspath(gate, full_path) || throw(ArgumentError("Gate does not have path $full_path"))
                result["$name_string/$pulse_name"] = (full_path, pulse_value)
            end
        else
            result[name_string] = (path, value)
        end
    end
    return result
end

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

_operator_matrix(operator::AbstractMatrix) = Matrix(operator)
_operator_matrix(operator::AbstractQuantumObject) = Matrix(operator.data)
_operator_matrix(operator) =
    throw(ArgumentError("Expected a matrix or quantum operator, received $(typeof(operator))"))

"""
    unitary_infidelity(target, actual)

Return `clamp(1 - abs2(tr(target' * actual))/d^2, 0, 1)` for equal square matrices
or quantum operators. The target must be unitary (tolerance `1e-8`). This is a
phase-insensitive process infidelity for unitary evolution; `actual` is not
independently checked for unitarity. It is not a state-transfer fidelity or a
noise-channel metric. For two unitaries, average gate infidelity is `d/(d+1)`
times this result.
"""
function unitary_infidelity(target, actual)
    target_matrix = _operator_matrix(target)
    actual_matrix = _operator_matrix(actual)
    size(target_matrix) == size(actual_matrix) || throw(DimensionMismatch(
        "Target size $(size(target_matrix)) does not match actual size $(size(actual_matrix))",
    ))
    size(target_matrix, 1) == size(target_matrix, 2) ||
        throw(DimensionMismatch("Unitary matrices must be square"))
    dimension = size(target_matrix, 1)
    identity = Matrix{eltype(target_matrix)}(I, dimension, dimension)
    isapprox(target_matrix' * target_matrix, identity; atol = 1e-8, rtol = 1e-8) ||
        throw(ArgumentError("Calibration target must be unitary"))
    return Float64(clamp(
        real(1 - abs2(tr(target_matrix' * actual_matrix)) / dimension^2),
        0,
        1,
    ))
end

"""
    gate_unitary(model::DeviceModel, gate::DeviceGate; kwargs...)

Evolve the identity over `[0, gate.duration]` with QuantumToolbox `sesolve` and
return the final retained-space propagator. Duration must be positive. Solver
kwargs pass through, with `progress_bar=false` by default; unsuccessful retcodes
raise an error. Units must obey the package's ℏ=1 Hamiltonian convention.
"""
function gate_unitary(model::DeviceModel, gate::DeviceGate; kwargs...)
    gate.duration > 0 || throw(ArgumentError("Gate evolution requires a positive duration"))
    options = merge((; progress_bar = false), (; kwargs...))
    solution = sesolve(
        numerical(model, gate),
        qeye_like(model.H),
        [zero(gate.duration), gate.duration];
        options...,
    )
    SciMLBase.successful_retcode(solution.retcode) ||
        error("Gate evolution failed with return code $(solution.retcode)")
    return solution.states[end]
end

"""
    gate_infidelity(model::DeviceModel, gate::DeviceGate, target)

Evaluate [`unitary_infidelity`](@ref) between `target` and
[`gate_unitary`](@ref). The target acts on the full retained model space.
This is the default calibration objective, not a computational-subspace
leakage metric or open-system fidelity.
"""
gate_infidelity(model::DeviceModel, gate::DeviceGate, target) =
    unitary_infidelity(target, gate_unitary(model, gate))

"""
    calibration_problem(model, gate, names, target; objective=gate_infidelity, kwargs...)

Build a standard `SciMLBase.OptimizationProblem` and an immutable gate rebuilder.
Select parameters using names from `parameters(gate)`, e.g. `["drive/amplitude"]`.
Additional keyword arguments (including `lb` and `ub`) are forwarded to
`OptimizationProblem`. The objective has signature `(model, candidate_gate, target)`
and must return a scalar loss. The default compares full retained-space unitaries.

Selected names must be a nonempty collection of unique discovery names with finite
real values. Returns a [`SciMLCalibrationSetup`](@ref) without modifying the input
gate. Algorithm selection and solver kwargs belong to [`calibrate`](@ref).
"""
function calibration_problem(model::DeviceModel, gate::DeviceGate, names, target;
                             objective = gate_infidelity, kwargs...)
    names isa AbstractString && throw(ArgumentError("Pass a list of parameter names, such as [\"drive/amplitude\"]"))
    names = collect(names)
    isempty(names) && throw(ArgumentError("At least one calibration parameter is required"))
    allunique(names) || throw(ArgumentError("Calibration parameter names must be unique"))
    initial_values, paths = _selected_calibration_parameters(gate, names)
    rebuild = values -> _with_calibration_values(gate, paths, values)
    loss(values, _) = objective(model, rebuild(values), target)
    problem = SciMLBase.OptimizationProblem(loss, initial_values; kwargs...)
    return SciMLCalibrationSetup(problem, rebuild)
end

"""
    calibrated_gate(setup::SciMLCalibrationSetup, values)

Reconstruct a gate from candidate parameter values in the selected order without
modifying the original gate. The vector length must match the calibration setup.
"""
calibrated_gate(setup::SciMLCalibrationSetup, values) = setup.rebuild(values)

"""
    calibrate(setup, algorithm; kwargs...)

Solve a calibration through the standard SciML interface and rebuild its gate.
All keyword arguments are forwarded directly to `SciMLBase.solve`.
Returns `(solution, gate)` with the candidate reconstructed from `solution.u`.
Inspect the solution retcode and objective before accepting the result; this wrapper
does not independently enforce optimization convergence. The original gate is
unchanged. For direct solve access use `setup.problem` and [`calibrated_gate`](@ref).
"""
function calibrate(setup::SciMLCalibrationSetup, algorithm; kwargs...)
    solution = SciMLBase.solve(setup.problem, algorithm; kwargs...)
    return solution, calibrated_gate(setup, solution.u)
end
