function expansive_tracking(A, B)
    nA = length(A)
    nB = length(B)
    nA <= nB || throw(ArgumentError("Cannot track more states than are available"))

    overlaps = [
        abs2(dot(A[i], B[j]))
        for i in 1:nA, j in 1:nB
    ]

    indices = Vector{Int}(undef, nA)
    overlap_metric = Vector{Float64}(undef, nA)

    for _ in 1:nA
        i, j = Tuple(argmax(overlaps))

        indices[i] = j
        overlap_metric[i] = overlaps[i, j]

        overlaps[i, :] .= -Inf
        overlaps[:, j] .= -Inf
    end

    return indices, overlap_metric
end