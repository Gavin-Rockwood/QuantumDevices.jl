"""A named builtin operator factory; configuration is portable without closure capture."""
struct BuiltinGenerator
    kind::Symbol
end
function (generator::BuiltinGenerator)(dimension; kwargs...)
    kind = generator.kind
    if kind in (:pauli_x,:pauli_y,:pauli_z)
        dimension == 2 || throw(ArgumentError("Pauli basis needs dimension 2"))
        return kind === :pauli_x ? sigmax() : kind === :pauli_y ? sigmay() : sigmaz()
    end
    kind === :destroy && return destroy(dimension)
    kind === :number && return num(dimension)
    kind === :charge && return num(charge_basis_dimension(dimension)) - (dimension-1) ÷ 2
    kind === :jump && return tunneling(charge_basis_dimension(dimension),1)
    throw(ArgumentError("Unknown builtin generator: $kind"))
end
