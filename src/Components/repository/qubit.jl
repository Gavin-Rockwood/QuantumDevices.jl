struct QubitSpec <: AbstractComponentSpec
    frequency::DeviceParameter
    dimension::Dimension
    operators
    hamiltonian
end

function QubitSpec(val;
    domain = realdomain(),
    fixed_frequency = false,
    description = "Qubit Frequency",
    metadata = Dict{Symbol, Any}(),
    hamiltonian = missing,
    extra_terms = nothing
    )
    frequency_parameter = DeviceParameter(
        ParamPath(:frequency),
        domain,
        fixed_frequency,
        true,
        val,
        description,
        _metadata_dict(metadata)
    )

    dimension = Dimension(2)
    operators = (
        identity = (; kwargs...) -> qeye(2),
        x = (; kwargs...) -> sigmax(),
        y = (; kwargs...) -> sigmay(),
        z = (; kwargs...) -> sigmaz(),
        p = (; kwargs...) -> sigmap(),
        m = (; kwargs...) -> sigmam(),
    )
    default_hamiltonian = param(frequency_parameter) * op(:z) / 2
    resolved_hamiltonian = hamiltonian === missing ? default_hamiltonian : hamiltonian
    resolved_hamiltonian = _component_hamiltonian(resolved_hamiltonian, extra_terms)
    return QubitSpec(frequency_parameter, dimension, operators, resolved_hamiltonian)
end
