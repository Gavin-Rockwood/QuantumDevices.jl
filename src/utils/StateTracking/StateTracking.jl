"""Result of tracking states through a sequence of candidate state collections."""
struct StateTrackingResult{S<:AbstractMatrix,O<:AbstractDict,D}
    states::S
    overlaps::Matrix{Float64}
    assignments::Matrix{Int}
    other_sorts::O
    diagnostics::D
end

"""Choose the largest-overlap candidate independently for every tracked state."""
greedy_assignment(overlaps::AbstractMatrix) =
    [argmax(view(overlaps, state, :)) for state in axes(overlaps, 1)]

"""
    track_states(state_history, states_to_track=first(state_history); kwargs...)

Track states by overlap with their selection at the preceding step. `method`
receives the complete overlap matrix for a step and returns one candidate index
per tracked state. Auxiliary candidate data follows the same step/candidate
layout as `state_history` through `other_sorts`.
"""
function track_states(
        state_history,
        states_to_track = first(state_history);
        other_sorts = Dict(),
        method = greedy_assignment,
        diagnostics::Symbol = :none,
        minimum_overlap::Real = 0.5)
    isempty(state_history) && error("State history cannot be empty.")
    diagnostics in (:none, :events, :full) ||
        error("Unknown diagnostics mode $diagnostics.")
    0 <= minimum_overlap <= 1 ||
        error("Minimum overlap must be between zero and one.")

    tracked = collect(_tracking_states(states_to_track))
    isempty(tracked) && error("At least one state must be tracked.")
    steps = length(state_history)
    states = Matrix{Any}(undef, length(tracked), steps)
    overlaps = Matrix{Float64}(undef, length(tracked), steps)
    assignments = Matrix{Int}(undef, length(tracked), steps)
    sorted = Dict(key => Matrix{Any}(undef, length(tracked), steps)
        for key in keys(other_sorts))
    collisions = NamedTuple[]
    low_overlaps = NamedTuple[]
    overlap_matrices = Matrix{Float64}[]

    previous = copy(tracked)
    for step in eachindex(state_history)
        candidates = collect(_tracking_states(state_history[step]))
        isempty(candidates) && error("State-history entry $step is empty.")
        overlap_matrix = Float64[
            abs2(dot(_state_vector(reference), _state_vector(candidate)))
            for reference in previous, candidate in candidates
        ]
        selected = Int.(collect(method(overlap_matrix)))
        length(selected) == length(tracked) ||
            error("Tracking method returned $(length(selected)) assignments for $(length(tracked)) states.")
        all(index -> checkbounds(Bool, candidates, index), selected) ||
            error("Tracking method returned an invalid candidate index.")

        assignments[:, step] = selected
        overlaps[:, step] = overlap_matrix[CartesianIndex.(axes(overlap_matrix, 1), selected)]
        states[:, step] = candidates[selected]
        previous = states[:, step]
        for key in keys(other_sorts)
            values = other_sorts[key][step]
            length(values) >= maximum(selected) ||
                error("Other sort $key has too few values at step $step.")
            sorted[key][:, step] = values[selected]
        end

        diagnostics == :none && continue
        for candidate in unique(selected)
            owners = findall(==(candidate), selected)
            length(owners) > 1 && push!(collisions,
                (; step, candidate, tracked_states = owners))
        end
        for state in eachindex(selected)
            overlaps[state, step] < minimum_overlap && push!(low_overlaps,
                (; step, tracked_state = state, overlap = overlaps[state, step]))
        end
        diagnostics == :full && push!(overlap_matrices, overlap_matrix)
    end

    diagnostic_data = (; collisions, low_overlaps, overlap_matrices)
    return StateTrackingResult(states, overlaps, assignments, sorted, diagnostic_data)
end

"""Generate an adiabatic state history and track its states and energies."""
function get_dressed_states(
        uncoupled_hamiltonian,
        interaction_hamiltonian,
        states_to_track = _tracking_eigenstates(uncoupled_hamiltonian).vectors;
        steps::Integer = 20,
        schedule::Symbol = :power,
        exponent::Real = 3.0,
        kwargs...)
    steps >= 2 || error("Dressed-state tracking requires at least two steps.")
    schedule in (:linear, :power) || error("Unknown dressing schedule $schedule.")
    exponent > 0 || error("Dressing schedule exponent must be positive.")

    progress = collect(range(0.0, 1.0; length = steps))
    strengths = schedule == :linear ? progress : progress .^ exponent
    eigensystems = [_tracking_eigenstates(
        uncoupled_hamiltonian + strength * interaction_hamiltonian,
    ) for strength in strengths]
    return track_states(
        [eigensystem.vectors for eigensystem in eigensystems],
        states_to_track;
        other_sorts = Dict(:energy => [real.(system.values) for system in eigensystems]),
        kwargs...,
    )
end

_tracking_states(states::AbstractMatrix) = eachcol(states)
_tracking_states(states) = states
_state_vector(state::AbstractVector) = state
_state_vector(state) = vec(state.data)

function _tracking_eigenstates(hamiltonian::AbstractMatrix)
    eigensystem = eigen(Hermitian(Matrix(hamiltonian)))
    return (values = eigensystem.values, vectors = eigensystem.vectors)
end

_tracking_eigenstates(hamiltonian) = eigenstates(hamiltonian)
