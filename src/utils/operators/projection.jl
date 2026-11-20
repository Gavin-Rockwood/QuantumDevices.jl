# Embedding: retained space → full input model space. Internal indices are one-based.
function selection(dimension, indices)
    return sum(
        basis(dimension, j - 1) * basis(length(indices), k - 1)' for (k, j) in enumerate(indices)
    )
end

ktensor(operators) = reduce(tensor, operators)
