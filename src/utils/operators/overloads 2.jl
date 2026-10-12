truncate(H::QuantumObject, P::QuantumObject) = QuantumObject(P.data' * H.data * P.data)
function truncate(H::QuantumObjectEvolution, P::QuantumObject)
    # Rectangular basis changes need separate input/output dimensions.
    return QuantumObjectEvolution(MatrixOperator(P.data') * H.data * MatrixOperator(P.data))
end
