function diagonalize(hamiltonian::QuantumObject)
    # Basis rotations can introduce roundoff asymmetry in Hermitian operators.
    isapprox(hamiltonian, hamiltonian'; rtol = 1e-12, atol = 1e-12) ||
        throw(ArgumentError("Hamiltonian must be Hermitian (rtol=atol=1e-12)"))
    eigensystem = eigenstates((hamiltonian + hamiltonian') / 2)
    # Flatten the eigenbasis while retaining the input tensor dimensions.
    rotation = QuantumObject(
        eigensystem.vectors;
        dims = (Tuple(hamiltonian.dims[1]), (length(eigensystem.values),)),
    )
    return Float64.(real.(eigensystem.values)), rotation
end

spectrum_operator(energies) = QuantumObject(Diagonal(energies))
