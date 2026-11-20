# Matcher contract: candidate indices and matched overlap amplitudes, both in
# reference order. track_states validates the nonempty input collections and
# candidate counts; amplitudes retain phase and any projected-reference norm loss.
# Assign in reference order, evaluating only candidates that remain available.
# Strict improvement preserves the first candidate on ties.
function quick_overlap_assignment(references, candidates)
    available = trues(length(candidates))
    assignment = Vector{Int}(undef, length(references))
    amplitudes = Vector{typeof(dot(first(references), first(candidates)))}(undef, length(references))
    for (i, reference) in enumerate(references)
        best = -Inf
        for j in eachindex(candidates)
            available[j] || continue
            amplitude = dot(reference, candidates[j])
            weight = abs2(amplitude)
            isfinite(weight) || throw(ArgumentError("nonfinite state overlaps"))
            if weight > best
                best = weight
                assignment[i] = j
                amplitudes[i] = amplitude
            end
        end
        available[assignment[i]] = false
    end
    return assignment, amplitudes
end

# Cached overlaps have rows already restored to tracked-label order.
function quick_overlap_assignment(overlaps::AbstractMatrix)
    available = trues(size(overlaps, 2))
    assignment = Vector{Int}(undef, size(overlaps, 1))
    amplitudes = Vector{eltype(overlaps)}(undef, size(overlaps, 1))
    for i in axes(overlaps, 1)
        best = -Inf
        for j in axes(overlaps, 2)
            available[j] || continue
            amplitude = overlaps[i, j]
            weight = abs2(amplitude)
            isfinite(weight) || throw(ArgumentError("nonfinite state overlaps"))
            if weight > best
                best = weight
                assignment[i] = j
                amplitudes[i] = amplitude
            end
        end
        available[assignment[i]] = false
    end
    return assignment, amplitudes
end
