"""Shared owner for tracked arrays, saved samples, keys, and full-path diagnostics."""
struct TrackingResult{K, S, O, I, A, D, L}
    keys::K
    saved_steps::Tuple{Vararg{Int}}
    total_steps::Int
    _states::S
    _other::O
    _indices::I
    _accuracy::A
    _diagnostics::D
    _lookup::L
    tracked::Bool
end

function TrackingResult(
    keys,
    saved_steps,
    total_steps,
    states,
    other::NamedTuple,
    indices,
    accuracy,
    diagnostics;
    tracked = true,
)
    keys = Tuple(keys)
    saved_steps = Tuple(saved_steps)
    rows = size(states, 1)
    samples = length(saved_steps)
    if !(
        all(i -> i isa Integer && !(i isa Bool), saved_steps) &&
        samples > 0 &&
        total_steps >= last(saved_steps) &&
        first(saved_steps) >= 1 &&
        all(i -> saved_steps[i - 1] < saved_steps[i], 2:samples)
    )
        throw(ArgumentError(
            "saved sample indices must be nonempty, increasing, and in the processed path",
        ))
    end

    tracked &&
        isempty(keys) &&
        throw(ArgumentError("only explicitly untracked results may have empty keys"))
    lookup = Dict(key => i for (i, key) in enumerate(keys))
    length(lookup) == length(keys) || throw(ArgumentError("tracking keys must be unique"))
    (tracked ? length(keys) == rows : isempty(keys)) || throw(
        DimensionMismatch(
            "tracking keys must match rows, or be empty for an explicitly untracked result",
        ),
    )
    expected = samples == 1 ? (rows,) : (rows, samples)
    size(states) == expected && size(indices) == expected ||
        throw(DimensionMismatch("saved states and indices must match rows and saved steps"))
    all(a -> a isa AbstractArray && size(a) == expected, values(other)) ||
        throw(DimensionMismatch("auxiliary data must match saved state dimensions"))
    if tracked
        accuracy isa NamedTuple &&
        all(a -> a isa AbstractArray && size(a) == expected, values(accuracy)) ||
            throw(DimensionMismatch("accuracy arrays must match saved state dimensions"))
        diagnostics isa AbstractVector && length(diagnostics) == rows ||
            throw(DimensionMismatch("final diagnostics must match rows"))
    else
        isnothing(accuracy) && isnothing(diagnostics) ||
            throw(ArgumentError("untracked results have no accuracy or diagnostics"))
    end
    return TrackingResult(
        keys,
        saved_steps,
        Int(total_steps),
        states,
        other,
        indices,
        accuracy,
        diagnostics,
        lookup,
        tracked,
    )
end

function tracking_row(result::TrackingResult, key)
    result.tracked || throw(
        ArgumentError(
            "occupation lookup is unavailable for this untracked eigensystem; enable recipe dressing",
        ),
    )
    return getfield(result, :_lookup)[key]
end
