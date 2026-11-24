# Reconstruct a transient reference map from the actual stored eigenstates.
function input_basis_map(model)
    dimension = length(model.states)
    return reduce(+, (ket * QT.basis(dimension, index - 1)' for (index, ket) in enumerate(model.states)))
end

# Compiler fixtures deliberately decouple the supplied spectrum from their test
# expressions. Build real models through build(); production has no fixture shim.
function fixture_numerical_model(recipe,spectrum,operators,hamiltonian,owners,eigenvectors,products)
    identity = QT.qeye(length(spectrum);dims=Tuple(eigenvectors.dims[1]))
    H = QuantumDevices.compile_hamiltonian(QuantumDevices.physical_terms(hamiltonian),products,
        recipe.resolved_parameters,Set{Symbolics.Num}(),identity)
    eigensystem = QuantumDevices.untracked_eigensystem(spectrum,eigenvectors,Tuple(sort!(collect(owners))))
    return NumericalModel(recipe,H,eigensystem,Dict{OperatorRef,QT.QuantumObject}(operators),hamiltonian,owners,products)
end
