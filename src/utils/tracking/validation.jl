function validate_tracking_threads(count; inherit=false)
    inherit && ismissing(count) && return count
    count isa Int && !(count isa Bool) && count > 0 ||
        throw(ArgumentError("tracking_threads must be a positive Int" * (inherit ? " or missing" : "")))
    return count
end

function tracking_method(method)
    (method === :quick || method === :global) ||
        throw(ArgumentError("method must be :quick or :global"))
    return method
end

function tracking_tolerances(atol, rtol, threshold)
    all(x -> x isa Real && isfinite(x) && x >= 0, (atol, rtol)) ||
        throw(ArgumentError("degeneracy tolerances must be finite nonnegative real values"))
    return threshold isa Real && isfinite(threshold) && 0 <= threshold <= 1 ||
           throw(ArgumentError("overlap_threshold must be finite and in [0,1]"))
end

function validate_kets(collection)
    isempty(collection) && throw(ArgumentError("state collections must be nonempty"))
    for ket in collection
        ket isa QuantumObject{Ket} ||
            throw(ArgumentError("state histories must contain QuantumObject kets"))
        isapprox(norm(ket), 1; atol = 1e-10, rtol = 1e-10) ||
            throw(ArgumentError("candidate states must be normalized"))
        ket.dimensions == first(collection).dimensions ||
            throw(DimensionMismatch("candidate ket dimensions differ"))
    end
    return nothing
end

function saved_sample_indices(selection, steps)
    selection === :all && return Tuple(1:steps)
    selection === :end && return (steps,)
    selection isa Union{AbstractVector, AbstractRange} &&
    !isempty(selection) &&
    all(i -> i isa Integer && !(i isa Bool) && 1 <= i <= steps, selection) &&
    issorted(selection) &&
    allunique(selection) || throw(
        ArgumentError(
            "save_steps must be :all, :end, or nonempty strictly increasing one-based sample indices",
        ),
    )
    return Tuple(Int.(selection))
end

function validate_tracking_transport(transport, previous, current)
    transport isa QuantumObject{Operator} ||
        throw(ArgumentError("transports must be QuantumObject linear maps"))
    transport.dims[2] == first(previous).dims[1] &&
    transport.dims[1] == first(current).dims[1] || throw(
        DimensionMismatch("transport must map previous ket coordinates into current ket coordinates"))
    return nothing
end

function prepare_tracking_inputs(state_history, other, transports; reference_count=nothing)
    history = [collect(collection) for collection in state_history]
    isempty(history) && throw(ArgumentError("state history must be nonempty"))
    foreach(validate_kets, history)
    n = isnothing(reference_count) ? length(first(history)) : reference_count
    steps = length(history)
    all(c -> length(c) >= n, history) || throw(
        ArgumentError("candidate collections cannot have fewer states than the tracked reference count"),
    )
    if !isnothing(transports)
        length(transports) == steps - 1 ||
            throw(DimensionMismatch("supply one transport per transition"))
        for t = 2:steps
            validate_tracking_transport(transports[t - 1], history[t - 1], history[t])
        end
    else
        all(c -> first(c).dimensions == first(first(history)).dimensions, history) ||
            throw(DimensionMismatch("different state bases require explicit transports"))
    end
    other isa NamedTuple ||
        throw(ArgumentError("other must be a NamedTuple of candidate-indexed histories"))
    auxiliary = map(other) do supplied
        length(supplied) == steps ||
            throw(DimensionMismatch("auxiliary histories must match input sample count"))
        values = [collect(v) for v in supplied]
        all(t -> length(values[t]) == length(history[t]), 1:steps) ||
            throw(DimensionMismatch("auxiliary data must have one value per candidate"))
        return values
    end
    if haskey(auxiliary, :energies)
        all(values -> all(x -> x isa Real && isfinite(x), values), auxiliary.energies) ||
            throw(ArgumentError("energies must be finite real values"))
    end
    return history, auxiliary
end
