"""
    make_qubit(name, ν)
    make_qubit(; name, ν, dimension=2)

Construct a two-level component with `H = ν * Z / 2` and operators `x`, `y`, `z`,
`p`, and `m` from QuantumToolbox. The input `ν` is an energy divided by ℏ
(angular frequency), with units reciprocal to simulation time. The dimension
is fixed at two. This Pauli convention assigns `+ν/2` to the first basis vector;
do not infer ground-state labels from its array index.
"""
function make_qubit(name, ν)
    parameters = (ν=ν,)
    operators = (x=sigmax(), y=sigmay(), z=sigmaz(), p=sigmap(), m=sigmam())
    hamiltonian = 0.5*param(:ν)*op(:z)
    return Component(name, parameters, operators, hamiltonian, "qubit", 2)
end

function make_qubit(; name, ν, dimension=2)
    dimension == 2 || throw(ArgumentError("Qubit dimension is fixed at 2"))
    return make_qubit(name, ν)
end

COMPONENT_CONSTRUCTORS["qubit"] = make_qubit
