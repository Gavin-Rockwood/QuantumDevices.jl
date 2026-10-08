"""
    state_amplitudes(references::AbstractDict, solution::QuantumToolbox.TimeEvolutionSol)

Compute reference-state amplitudes from `solution.states` using the generic
history helper. The returned vectors correspond to `solution.times_states`,
not necessarily `solution.times` (the expectation-value sampling grid).
The saved states must be kets; density-matrix solutions are rejected.
"""
state_amplitudes(references::AbstractDict, solution::QuantumToolbox.TimeEvolutionSol) =
    state_amplitudes(references, solution.states)

"""
    sesolve(model::DeviceModel, gate::DeviceGate, initial, times;
            energy_shift=:auto, dense=false, kwargs...)

Evolve a ket, rectangular operator, or vector of kets using the model and gate.
Frequency units are converted with `numerical(...; scalar=2pi)`, and pulse event
stops are supplied automatically. `energy_shift=:auto` subtracts the mean idle
energy of the input states during integration and restores global phase at each
saved time. An explicit offset is in frequency units; `energy_shift=0` disables it.
`dense=true` constructs dense Hamiltonian matrices before integration; the default
preserves their usual storage.
Returns the native QuantumToolbox solution and forwards native solver options.
User integrator callbacks observe the centered state during integration.
"""
function QuantumToolbox.sesolve(model::DeviceModel, gate::DeviceGate,
    initial::Union{QuantumObject,AbstractVector{<:QuantumObject}}, times::AbstractVector;
    energy_shift=:auto, dense::Bool=false, kwargs...)
    isempty(times) && throw(ArgumentError("Evolution times cannot be empty"))
    columns = initial isa AbstractVector ? _state_columns(initial) : initial
    shift = _energy_shift(model.H, columns, energy_shift)
    evaluate = _prepared_gate_numerical(model, gate; scalar=2pi, energy_shift=shift, dense)
    stops = filter(t -> first(times) < t < last(times), pulse_tstops(gate))
    options = merge((; tstops=stops), (; kwargs...))
    solution = QuantumToolbox.sesolve(evaluate(gate), columns, times; options...)
    return _restore_global_phase!(solution, 2pi*shift, first(times))
end

# Native-Hamiltonian overloads retain angular-frequency units and native return types.
# Batched sesolve centers automatically. sesolveProblem remains a native constructor.
QuantumToolbox.sesolve(H::Union{AbstractQuantumObject{Operator},Tuple},
    states::AbstractVector{<:QuantumObject}, times::AbstractVector; energy_shift=:auto, kwargs...) =
    _sesolve_centered(H, _state_columns(states), times; energy_shift, kwargs...)

QuantumToolbox.sesolveProblem(H::Union{AbstractQuantumObject{Operator},Tuple},
    states::AbstractVector{<:QuantumObject}, times::AbstractVector; kwargs...) =
    QuantumToolbox.sesolveProblem(H, _state_columns(states), times; kwargs...)

function get_gate_matrix(solution::QuantumToolbox.TimeEvolutionSol, output_states;
    frame=:lab, idle_hamiltonian=nothing)
    SciMLBase.successful_retcode(solution.retcode) || error("Cannot project an unsuccessful evolution")
    isempty(solution.states) && throw(ArgumentError("Solution has no saved operator states"))
    columns = last(solution.states)
    isoper(columns) || throw(ArgumentError("Gate matrices require rectangular operator evolution"))
    output = _orthonormal_columns(output_states)
    output.dimensions.to == columns.dimensions.to ||
        throw(DimensionMismatch("Output basis and evolved states have different dimensions"))
    projection = _prepare_gate_projection(Matrix(output.data), frame, idle_hamiltonian)
    elapsed = last(solution.times_states) - first(solution.times)
    return _project_gate_columns(columns.data, projection, elapsed)
end
