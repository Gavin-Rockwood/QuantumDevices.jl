function projection(N::Integer, inds)
    N > 0 || throw(ArgumentError("Parent dimension must be positive"))
    indices = collect(inds)
    isempty(indices) && throw(ArgumentError("Projection must retain at least one state"))
    all(i -> i isa Integer && 1 <= i <= N, indices) || throw(ArgumentError("Invalid projection indices"))
    allunique(indices) || throw(ArgumentError("Projection indices must be unique"))
    return QuantumObject(sparse(indices, eachindex(indices), ones(length(indices)), N, length(indices)))
end

truncate(H::QuantumObject, P::QuantumObject) = QuantumObject(P.data' * H.data * P.data)
function truncate(H::QuantumObjectEvolution, P::QuantumObject)
    # Rectangular basis changes need separate input/output dimensions.
    return QuantumObjectEvolution(MatrixOperator(P.data') * H.data * MatrixOperator(P.data))
end
