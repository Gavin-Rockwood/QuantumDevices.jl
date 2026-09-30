function quick_tracking(A, B)
    nA = length(A)
    nB = length(B)
    nA <= nB || throw(ArgumentError("Cannot track more states than are available"))

    indices = Vector{Int}(undef, nA)
    overlap_metric = Vector{Float64}(undef, nA)
    available = collect(1:nB)

    for i in 1:nA
        overlaps = [abs2(dot(A[i], B[j])) for j in available]

        k = argmax(overlaps)
        j = available[k]

        indices[i] = j
        overlap_metric[i] = overlaps[k]

        deleteat!(available, k)
    end

    return indices, overlap_metric
end