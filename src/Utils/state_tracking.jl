function maximum_overlap_assignment(overlaps::AbstractMatrix{<:Real})
    tracked_count, candidate_count = size(overlaps)
    tracked_count <= candidate_count ||
        error("Cannot assign $tracked_count tracked states to $candidate_count candidate states.")
    tracked_count == 0 && return Int[]

    costs = -Float64.(overlaps)
    row_potential = zeros(Float64, tracked_count + 1)
    column_potential = zeros(Float64, candidate_count + 1)
    matched_row = zeros(Int, candidate_count + 1)
    previous_column = zeros(Int, candidate_count + 1)

    for row in 1:tracked_count
        matched_row[1] = row
        column = 1
        minimum_cost = fill(Inf, candidate_count + 1)
        used = falses(candidate_count + 1)

        while true
            used[column] = true
            current_row = matched_row[column]
            delta = Inf
            next_column = 1
            for candidate_column in 2:(candidate_count + 1)
                used[candidate_column] && continue
                reduced_cost =
                    costs[current_row, candidate_column - 1] -
                    row_potential[current_row + 1] -
                    column_potential[candidate_column]
                if reduced_cost < minimum_cost[candidate_column]
                    minimum_cost[candidate_column] = reduced_cost
                    previous_column[candidate_column] = column
                end
                if minimum_cost[candidate_column] < delta
                    delta = minimum_cost[candidate_column]
                    next_column = candidate_column
                end
            end

            for candidate_column in 1:(candidate_count + 1)
                if used[candidate_column]
                    row_potential[matched_row[candidate_column] + 1] += delta
                    column_potential[candidate_column] -= delta
                else
                    minimum_cost[candidate_column] -= delta
                end
            end
            column = next_column
            matched_row[column] == 0 && break
        end

        while true
            previous = previous_column[column]
            matched_row[column] = matched_row[previous]
            column = previous
            column == 1 && break
        end
    end

    assignment = zeros(Int, tracked_count)
    for column in 2:(candidate_count + 1)
        row = matched_row[column]
        row == 0 && continue
        assignment[row] = column - 1
    end
    return assignment
end

function track_state_history(
        state_history,
        reference_states::AbstractMatrix;
        histories = Dict(),
        minimum_overlap::Real = 0.5,
        on_low_overlap::Symbol = :warn)
    isempty(state_history) && error("State history cannot be empty.")
    0 <= minimum_overlap <= 1 ||
        error("Minimum overlap must be between zero and one.")
    on_low_overlap in (:warn, :error, :ignore) ||
        error("Unknown low-overlap policy $on_low_overlap.")

    tracked_count = size(reference_states, 2)
    steps = length(state_history)
    overlaps = Matrix{Float64}(undef, tracked_count, steps)
    assignments = Matrix{Int}(undef, tracked_count, steps)
    tracked_histories = _initialize_tracked_histories(histories, tracked_count, steps)
    previous_states = Matrix{ComplexF64}(reference_states)

    for step in 1:steps
        candidates = state_history[step]
        size(candidates, 1) == size(reference_states, 1) ||
            error("State-history entry $step has an incompatible state dimension.")
        overlap_matrix = abs2.(adjoint(previous_states) * candidates)
        assignment = maximum_overlap_assignment(overlap_matrix)
        selected_overlaps = [
            overlap_matrix[index, assignment[index]]
            for index in 1:tracked_count
        ]
        overlaps[:, step] = selected_overlaps
        assignments[:, step] = assignment
        _record_tracked_histories!(
            tracked_histories,
            histories,
            assignment,
            step,
        )
        previous_states = candidates[:, assignment]
        _handle_low_tracking_overlap(
            selected_overlaps,
            step,
            minimum_overlap,
            on_low_overlap,
        )
    end

    return (
        states = previous_states,
        overlaps = overlaps,
        assignments = assignments,
        histories = tracked_histories,
    )
end

function track_state_history(
        state_history::AbstractVector{<:AbstractVector},
        reference_states::AbstractDict;
        kwargs...)
    isempty(state_history) && error("State history cannot be empty.")
    if isempty(reference_states)
        generated_references = Dict(
            index => state
            for (index, state) in enumerate(first(state_history))
        )
        return track_state_history(
            state_history,
            generated_references;
            kwargs...,
        )
    end

    keys_to_track = collect(keys(reference_states))
    reference_matrix = hcat(
        (_state_vector(reference_states[key]) for key in keys_to_track)...,
    )
    matrix_history = [
        hcat((_state_vector(state) for state in slice)...)
        for slice in state_history
    ]
    tracking = track_state_history(matrix_history, reference_matrix; kwargs...)
    final_slice = last(state_history)
    final_assignment = tracking.assignments[:, end]

    return Dict(
        key => final_slice[final_assignment[index]]
        for (index, key) in enumerate(keys_to_track)
    )
end

function adiabatic_sweep(
        hamiltonian,
        reference_states::AbstractMatrix;
        lambdas = range(0.0, 1.0; length = 20),
        histories = Dict(),
        minimum_overlap::Real = 0.5,
        on_low_overlap::Symbol = :warn)
    lambda_values = collect(lambdas)
    isempty(lambda_values) && error("Adiabatic sweep requires at least one lambda value.")
    0 <= minimum_overlap <= 1 ||
        error("Minimum overlap must be between zero and one.")
    on_low_overlap in (:warn, :error, :ignore) ||
        error("Unknown low-overlap policy $on_low_overlap.")

    initial_hamiltonian = hamiltonian(first(lambda_values))
    size(reference_states, 1) == size(initial_hamiltonian, 1) ||
        error("Reference states and Hamiltonian dimensions do not match.")

    tracked_count = size(reference_states, 2)
    steps = length(lambda_values)
    overlaps = Matrix{Float64}(undef, tracked_count, steps)
    assignments = Matrix{Int}(undef, tracked_count, steps)
    energy_history = Matrix{Float64}(undef, tracked_count, steps)
    tracked_histories = _initialize_tracked_histories(histories, tracked_count, steps)
    previous_states = Matrix{ComplexF64}(reference_states)
    final_states = previous_states
    final_energies = zeros(Float64, tracked_count)

    for step in 1:steps
        current_hamiltonian =
            step == 1 ? initial_hamiltonian : hamiltonian(lambda_values[step])
        size(current_hamiltonian) == size(initial_hamiltonian) ||
            error("Hamiltonian at sweep step $step has an incompatible dimension.")
        eigensystem = _tracking_eigenstates(current_hamiltonian)
        overlap_matrix = abs2.(adjoint(previous_states) * eigensystem.vectors)
        assignment = maximum_overlap_assignment(overlap_matrix)
        selected_states = eigensystem.vectors[:, assignment]
        selected_overlaps = [
            overlap_matrix[index, assignment[index]]
            for index in 1:tracked_count
        ]

        overlaps[:, step] = selected_overlaps
        assignments[:, step] = assignment
        energy_history[:, step] = real.(eigensystem.values[assignment])
        _record_sweep_histories!(
            tracked_histories,
            histories,
            assignment,
            step,
            lambda_values[step],
            eigensystem,
        )
        previous_states = selected_states
        final_states = selected_states
        final_energies = energy_history[:, step]

        _handle_low_tracking_overlap(
            selected_overlaps,
            step,
            minimum_overlap,
            on_low_overlap,
        )
    end

    return (
        states = final_states,
        energies = final_energies,
        overlaps = overlaps,
        assignments = assignments,
        lambdas = lambda_values,
        energy_history = energy_history,
        histories = tracked_histories,
    )
end

function get_dressed_states(
        uncoupled_hamiltonian,
        interaction_hamiltonian;
        kwargs...)
    references = _tracking_eigenstates(uncoupled_hamiltonian).vectors
    return get_dressed_states(
        uncoupled_hamiltonian,
        interaction_hamiltonian,
        references;
        kwargs...,
    )
end

function get_dressed_states(
        uncoupled_hamiltonian,
        interaction_hamiltonian,
        reference_states::AbstractMatrix;
        steps::Integer = 20,
        schedule::Symbol = :power,
        exponent::Real = 3.0,
        histories = Dict(),
        minimum_overlap::Real = 0.5,
        on_low_overlap::Symbol = :warn)
    steps >= 2 || error("Dressed-state tracking requires at least two steps.")
    schedule in (:linear, :power) ||
        error("Unknown dressing schedule $schedule.")
    exponent > 0 || error("Dressing schedule exponent must be positive.")

    progress = collect(range(0.0, 1.0; length = steps))
    lambdas = schedule == :linear ? progress : progress .^ exponent

    return adiabatic_sweep(
        lambda -> uncoupled_hamiltonian + lambda * interaction_hamiltonian,
        reference_states;
        lambdas = lambdas,
        histories = histories,
        minimum_overlap = minimum_overlap,
        on_low_overlap = on_low_overlap,
    )
end

function _initialize_tracked_histories(histories, tracked_count, steps)
    return Dict(
        key => Matrix{Any}(undef, tracked_count, steps)
        for key in keys(histories)
    )
end

function _record_tracked_histories!(
        tracked_histories,
        histories,
        assignment,
        step)
    for (key, history) in pairs(histories)
        length(history) >= step ||
            error("History $key does not contain sweep step $step.")
        _record_candidate_values!(
            tracked_histories[key],
            history[step],
            assignment,
            step,
            key,
        )
    end
end

function _record_sweep_histories!(
        tracked_histories,
        histories,
        assignment,
        step,
        lambda,
        eigensystem)
    for (key, history) in pairs(histories)
        candidate_values = if applicable(history, lambda, eigensystem)
            history(lambda, eigensystem)
        else
            length(history) >= step ||
                error("History $key does not contain sweep step $step.")
            history[step]
        end
        _record_candidate_values!(
            tracked_histories[key],
            candidate_values,
            assignment,
            step,
            key,
        )
    end
end

function _record_candidate_values!(
        destination,
        candidate_values,
        assignment,
        step,
        key)
    length(candidate_values) >= maximum(assignment; init = 0) ||
        error("History $key has too few candidate values at sweep step $step.")
    destination[:, step] = candidate_values[assignment]
end

_tracking_matrix(hamiltonian::AbstractMatrix) = hamiltonian
_tracking_matrix(hamiltonian) = hamiltonian.data
_state_vector(state::AbstractVector) = state
_state_vector(state) = vec(state.data)

function _tracking_eigenstates(hamiltonian::AbstractMatrix)
    eigensystem = eigen(Hermitian(Matrix(hamiltonian)))
    return (values = eigensystem.values, vectors = eigensystem.vectors)
end

_tracking_eigenstates(hamiltonian) = eigenstates(hamiltonian)

function _handle_low_tracking_overlap(
        overlaps,
        step,
        minimum_overlap,
        on_low_overlap)
    low = findall(value -> value < minimum_overlap, overlaps)
    isempty(low) && return
    message = "Dressed-state tracking overlap below $minimum_overlap at step $step for tracked states $(collect(low))."
    if on_low_overlap == :warn
        @warn message
    elseif on_low_overlap == :error
        error(message)
    end
end
