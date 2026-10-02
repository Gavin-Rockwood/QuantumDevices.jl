_operator_matrix(operator::AbstractMatrix) = Matrix(operator)
_operator_matrix(operator::AbstractQuantumObject) = Matrix(operator.data)
_operator_matrix(operator) =
    throw(ArgumentError("Expected a matrix or quantum operator, received $(typeof(operator))"))

"""
    unitary_fidelity(target, actual)

Return `clamp(abs2(tr(target' * actual))/d^2, 0, 1)` for equal square matrices
or quantum operators. The target must be unitary (tolerance `1e-8`); `actual`
is not independently checked. This phase-insensitive process fidelity compares
full unitaries, not state transfers or dissipative channels.
"""
function unitary_fidelity(target, actual)
    target_matrix = _operator_matrix(target)
    actual_matrix = _operator_matrix(actual)
    size(target_matrix) == size(actual_matrix) || throw(DimensionMismatch(
        "Target size $(size(target_matrix)) does not match actual size $(size(actual_matrix))",
    ))
    size(target_matrix, 1) == size(target_matrix, 2) ||
        throw(DimensionMismatch("Unitary matrices must be square"))
    dimension = size(target_matrix, 1)
    identity = Matrix{eltype(target_matrix)}(I, dimension, dimension)
    isapprox(target_matrix' * target_matrix, identity; atol = 1e-8, rtol = 1e-8) ||
        throw(ArgumentError("Target must be unitary"))
    return Float64(clamp(
        real(abs2(tr(target_matrix' * actual_matrix)) / dimension^2),
        0,
        1,
    ))
end

"""
    unitary_infidelity(target, actual)

Return `1 - unitary_fidelity(target, actual)`. For two unitaries, average gate
infidelity is `d/(d+1)` times this process infidelity. The target must be unitary;
the actual operator is not independently checked for unitarity.
"""
unitary_infidelity(target, actual) = 1 - unitary_fidelity(target, actual)
