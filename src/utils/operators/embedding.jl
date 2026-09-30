function wrap(operator, i, dims)
    1 <= i <= length(dims) || throw(ArgumentError("Invalid component index"))
    size(operator) == (dims[i], dims[i]) || throw(DimensionMismatch("Operator and component dimensions differ"))
    return tensor((i == j ? operator : qeye(dims[j]) for j in eachindex(dims))...)
end
