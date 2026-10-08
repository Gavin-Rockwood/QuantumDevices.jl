_operator_matrix(operator::AbstractMatrix) = Matrix(operator)
_operator_matrix(operator::AbstractQuantumObject) = Matrix(operator.data)
_operator_matrix(operator) =
    throw(ArgumentError("Expected a matrix or quantum operator, received $(typeof(operator))"))

"""
    unitary_fidelity(target, actual; include_phases=true)

Return `clamp(abs2(tr(target' * actual))/d^2, 0, 1)` for equal square matrices
or quantum operators. The target must be unitary (tolerance `1e-8`); `actual`
is not independently checked and may be a projected, leaking gate matrix.
The default includes relative phases while ignoring a common global phase.
With `include_phases=false`, return the average squared probability overlap
`sum(sum(abs.(target[:,j]) .* abs.(actual[:,j]))^2 for j)/d` instead.
This probability score ignores all entry phases and retains leakage penalties;
it is not coherent process fidelity. Neither mode normalizes `actual`.
"""
function unitary_fidelity(target, actual; include_phases::Bool=true)
    target_matrix = _validated_unitary_target(target)
    return _validated_matrix_fidelity(target_matrix, _operator_matrix(actual), include_phases)
end

function _validated_unitary_target(target)
    target_matrix = _operator_matrix(target)
    size(target_matrix, 1) == size(target_matrix, 2) ||
        throw(DimensionMismatch("Unitary matrices must be square"))
    dimension = size(target_matrix, 1)
    dimension > 0 || throw(ArgumentError("Target must be nonempty"))
    identity = Matrix{eltype(target_matrix)}(I, dimension, dimension)
    isapprox(target_matrix' * target_matrix, identity; atol = 1e-8, rtol = 1e-8) ||
        throw(ArgumentError("Target must be unitary"))
    return target_matrix
end

function _validated_matrix_fidelity(target, actual, include_phases::Bool)
    size(target) == size(actual) || throw(DimensionMismatch(
        "Target size $(size(target)) does not match actual size $(size(actual))"))
    dimension = size(target, 1)
    score = if include_phases
        abs2(dot(target, actual)) / dimension^2
    else
        sum(axes(target, 2)) do j
            abs2(sum(i -> abs(target[i, j]) * abs(actual[i, j]), axes(target, 1)))
        end / dimension
    end
    return Float64(clamp(real(score), 0, 1))
end

"""
    unitary_infidelity(target, actual; include_phases=true)

Return `1 - unitary_fidelity(target, actual)`. For two unitaries, average gate
infidelity is `d/(d+1)` times this process infidelity when phases are included.
With `include_phases=false`, return one minus the probability-overlap score.
The target must be unitary; the actual operator is not checked for unitarity.
"""
unitary_infidelity(target, actual; include_phases::Bool=true) =
    1 - unitary_fidelity(target, actual; include_phases)
