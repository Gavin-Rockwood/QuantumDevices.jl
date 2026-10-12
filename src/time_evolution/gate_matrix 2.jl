function _state_columns(states::AbstractVector{<:QuantumObject})
    isempty(states) && throw(ArgumentError("Pass at least one initial ket"))
    all(isket, states) || throw(ArgumentError("Initial states must be kets"))
    dimensions = first(states).dimensions
    all(s -> s.dimensions == dimensions, states) ||
        throw(DimensionMismatch("All initial kets must have the same quantum dimensions"))
    data = hcat((state.data for state in states)...)
    return QuantumObject(data; dims=QuantumToolbox.Dimensions(
        dimensions.to, QuantumToolbox.Space(length(states))))
end

function _orthonormal_columns(states)
    states isa AbstractVector{<:QuantumObject} ||
        throw(ArgumentError("Pass an ordered vector of QuantumToolbox kets"))
    columns = _state_columns(states)
    matrix = Matrix(columns.data)
    isapprox(matrix' * matrix, Matrix{eltype(matrix)}(I, length(states), length(states));
        atol=1e-8, rtol=1e-8) || throw(ArgumentError("Gate basis states must be orthonormal"))
    return columns
end

function _gate_bases(states, output_states)
    initial = _orthonormal_columns(states)
    output = _orthonormal_columns(output_states)
    size(initial) == size(output) || throw(DimensionMismatch(
        "Input and output gate bases must have equal dimensions and state counts"))
    initial.dimensions.to == output.dimensions.to ||
        throw(DimensionMismatch("Input and output bases have different quantum dimensions"))
    return initial, Matrix(output.data)
end

function _prepare_gate_projection(output, frame, idle_hamiltonian)
    frame in (:lab, :interaction) || throw(ArgumentError("Frame must be :lab or :interaction"))
    frame === :lab && return (; left=output', energies=nothing, right=nothing)
    idle_hamiltonian === nothing && throw(ArgumentError(
        "Interaction-frame comparison requires idle_hamiltonian in frequency units"))
    idle = _operator_matrix(idle_hamiltonian)
    size(idle) == (size(output, 1), size(output, 1)) ||
        throw(DimensionMismatch("Idle Hamiltonian and output basis dimensions differ"))
    ishermitian(idle) || throw(ArgumentError("Idle Hamiltonian must be Hermitian"))
    system = eigen(Hermitian(idle))
    return (; left=output' * system.vectors, energies=system.values, right=system.vectors')
end

function _project_gate_columns(columns, projection, duration)
    projection.energies === nothing && return projection.left * columns
    phases = cis.(2pi * duration .* projection.energies)
    return projection.left * (phases .* (projection.right * columns))
end

# H is already in angular-frequency units; numerical folds 2π into its matrices.
function _evolve_gate_columns(H, initial, duration; phase_shift=0, kwargs...)
    isfinite(duration) && duration > 0 || throw(ArgumentError("Evolution duration must be finite and positive"))
    options = merge((; progress_bar=false, saveat=[duration], save_start=false), (; kwargs...))
    result = sesolve(H, initial, [zero(duration), duration]; options...)
    SciMLBase.successful_retcode(result.retcode) || error("Gate evolution failed with retcode $(result.retcode)")
    isempty(result.states) && throw(ArgumentError("Gate evolution must save its final state"))
    last(result.times_states) == duration || throw(ArgumentError("Gate evolution must save at the gate duration"))
    _restore_global_phase!(result, phase_shift, zero(duration))
    return last(result.states).data
end

"""
    get_gate_matrix(model, gate, states; output_states=states, frame=:lab, kwargs...)
    get_gate_matrix(solution, output_states; frame=:lab, idle_hamiltonian=nothing)

Return `C' * Ψ(T)`, projecting simultaneously evolved input-state columns onto
the ordered output basis `C`. Bases are ordered vectors of orthonormal
QuantumToolbox kets. The model/gate form evolves with explicit `2π * H` and
forwards solver keywords. Its output basis defaults to the input basis.
The model/gate form is a compatibility alias for `numerical(model, gate, states; kwargs...)`
and returns a `QuantumObject` operator with dimensions matching the selected basis.

The solution form returns a matrix by projecting the last saved rectangular operator; its columns
retain the initial-state order. `frame=:lab` preserves lab-frame phases.
`frame=:interaction` removes idle evolution before projection, using the model's
idle Hamiltonian or the explicitly supplied `idle_hamiltonian` in frequency
units. Elapsed time is measured from `first(solution.times)`. Frame conversion
does not change integration. Projected columns are never renormalized, so
leakage remains visible. `1 - real(sum(abs2, M))/size(M, 2)` is mean leakage
out of the chosen output subspace for normalized evolved columns.
"""
function get_gate_matrix end
