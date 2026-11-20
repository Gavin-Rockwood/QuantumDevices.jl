"""An ephemeral construction-basis capture of a persistent recipe.

H is the stored reference QuantumObject, assembled and cleaned once at build
using the default numerical cleanup threshold. Operators, complete products, and
energy-ordered eigenstates share its coordinates. Explicit downstream projection
creates a new basis. Arrays retain ordinary mutability; H is an independent
captured reference, not a live view recomputed after cache mutations.
"""
struct NumericalModel{D <: DressedStates}
    recipe::ResolvedRecipe
    H::QuantumObject
    eigensystem::D
    operators::Dict{OperatorRef,QuantumObject}
    hamiltonian::Sym
    owners::Set{Symbol}
    operator_products::Dict{Tuple,QuantumObject}
end
