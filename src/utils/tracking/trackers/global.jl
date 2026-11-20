# Matcher contract: candidate indices and matched overlap amplitudes, both in
# reference order. track_states validates the nonempty input collections and
# candidate counts; amplitudes retain phase and any projected-reference norm loss.
# Return amplitudes as well as indices so the tracking loop can align phases.
function global_overlap_assignment(references, candidates)
    amplitudes = [dot(reference, candidate) for reference in references, candidate in candidates]
    assignment = maximum_overlap_assignment(abs2.(amplitudes))
    return assignment, [amplitudes[i, j] for (i, j) in enumerate(assignment)]
end

# Shortest augmenting paths for a rectangular maximum-weight assignment.
# Ascending scans and strict improvements give deterministic index tie-breaking.
function maximum_overlap_assignment(weights::AbstractMatrix{<:Real})
    n, m = size(weights)
    n > 0 && m >= n || throw(ArgumentError("assignment needs at least one row and no fewer columns"))
    all(isfinite, weights) || throw(ArgumentError("assignment weights must be finite"))
    u = zeros(n + 1)
    v = zeros(m + 1)
    matched = zeros(Int, m + 1)
    way = zeros(Int, m + 1)
    for row = 1:n
        matched[1] = row
        column = 1
        distances = fill(Inf, m + 1)
        used = falses(m + 1)
        while true
            used[column] = true
            i = matched[column]
            delta = Inf
            next = 0
            for j = 2:(m + 1)
                used[j] && continue
                cost = -Float64(weights[i, j - 1]) - u[i + 1] - v[j]
                if cost < distances[j]
                    distances[j] = cost
                    way[j] = column
                end
                if distances[j] < delta
                    delta = distances[j]
                    next = j
                end
            end
            for j = 1:(m + 1)
                if used[j]
                    u[matched[j] + 1] += delta
                    v[j] -= delta
                else
                    distances[j] -= delta
                end
            end
            column = next
            matched[column] == 0 && break
        end
        while column != 1
            previous = way[column]
            matched[column] = matched[previous]
            column = previous
        end
    end
    assignment = zeros(Int, n)
    for column = 2:(m + 1)
        matched[column] > 0 && (assignment[matched[column]] = column - 1)
    end
    return assignment
end

function global_overlap_assignment(amplitudes::AbstractMatrix)
    assignment = maximum_overlap_assignment(abs2.(amplitudes))
    return assignment, [amplitudes[i, j] for (i, j) in enumerate(assignment)]
end
