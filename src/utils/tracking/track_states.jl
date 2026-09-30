"""
    QuantumDevices.TrackingResult(states, others, confidence)

Result of [`track_states`](@ref) or [`get_dressed_states`](@ref).
For `track_states`, `states` maps labels to vectors of states at returned steps,
`others` holds correspondingly reordered auxiliary histories, and `confidence`
holds cumulative products of matched squared overlaps. For `get_dressed_states`,
these are final label-to-state, label-to-energy, and label-to-confidence dictionaries.
Confidence is a continuity diagnostic, not a calibrated probability of correctness.
"""
struct TrackingResult
    states
    others
    confidence
end

"""
    track_states(state_history, states_to_track=nothing;
                 other_sorts=Dict(), method=QuantumDevices.expansive_tracking,
                 return_steps=:all)

Match states through a nonempty sequence of eigenstate collections using squared
overlaps. Initial states may be a dictionary of labels or an indexed collection;
by default the first history entry supplies the labels. Each match is one-to-one.

`other_sorts` maps names to histories indexed like `state_history` (for example,
energies). `return_steps` is `:all`, `:end`, or valid step indices. Tracking still
passes through every intermediate step. A custom `method(A, B)` must return
matched indices into `B` and squared overlaps, in `A` order.

Returns `QuantumDevices.TrackingResult`. Both built-in methods are greedy overlap
assignments; neither promises a globally optimal assignment or a unique choice
inside a degenerate subspace.
"""
function track_states(
    state_history,
    states_to_track = nothing;
    other_sorts = Dict(),
    method = expansive_tracking,
    return_steps = :all,
)
    isempty(state_history) && throw(ArgumentError("State history cannot be empty"))
    states_to_track === nothing && (states_to_track = first(state_history))
    states_to_track_keys = states_to_track isa AbstractDict ?
        collect(keys(states_to_track)) :
        collect(eachindex(states_to_track))

    if states_to_track isa AbstractDict
        states_to_track = [states_to_track[key] for key in states_to_track_keys]
    end

    if return_steps === :all
        return_steps = eachindex(state_history)
    elseif return_steps === :end
        return_steps = [lastindex(state_history)]
    end
    return_steps isa Symbol && throw(ArgumentError("Unknown return_steps option"))
    all(i -> i in eachindex(state_history), return_steps) || throw(ArgumentError("Invalid return step"))
    for history in values(other_sorts)
        length(history) == length(state_history) || throw(DimensionMismatch("Auxiliary history length differs"))
        all(length(a) == length(b) for (a, b) in zip(history, state_history)) || throw(DimensionMismatch("Auxiliary state counts differ"))
    end

    indices = []
    overlap_metrics = []
    logconfidence = []

    # Track through the entire history
    for i in eachindex(state_history)
        if i == firstindex(state_history)
            idxs, overlap_metric = method(
                states_to_track,
                state_history[i],
            )

            logconf = log.(overlap_metric)
        else
            idxs, overlap_metric = method(
                state_history[i-1][indices[i-1]],
                state_history[i],
            )

            logconf = logconfidence[end] .+ log.(overlap_metric)
        end

        push!(indices, idxs)
        push!(overlap_metrics, overlap_metric)
        push!(logconfidence, logconf)
    end

    # Return only requested state-history steps
    states = Dict()

    for (k, key) in enumerate(states_to_track_keys)
        states[key] = [
            state_history[i][indices[i][k]]
            for i in return_steps
        ]
    end

    # Propagate any other state-indexed quantities
    others = Dict()

    for key in keys(other_sorts)
        others[key] = [
            other_sorts[key][i][indices[i]]
            for i in return_steps
        ]
    end

    # Cumulative confidence at each returned step
    confidence = [
        exp.(logconfidence[i])
        for i in return_steps
    ]

    return TrackingResult(states, others, confidence)
end
