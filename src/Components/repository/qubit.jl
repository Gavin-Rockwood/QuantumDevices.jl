struct QubitSpec <: AbstractComponentSpec
    frequency::DeviceParameter
    dimension::DeviceParameter
    operators
    hamiltonian
end

function QubitSpec(val;
    domain = realdomain(),
    fixed_frequency = false,
    description = "Qubit Frequency",
    metadata = Dict{Symbol, Any}(),
    hamiltonian = param(:frequency)*op(:z) / 2,
    extra_terms = nothing
    )
    frequency = DeviceParameter(
        ParamPath(:frequency),
        domain,
        fixed_frequency,
        true,
        val,
        description,
        _metadata_dict(metadata)
    )

    dimension = DeviceParameter(
        ParamPath(:dimension),
        integerrange(2,2),
        false,
        true,
        2,
        "Qubit Dimension",
        _metadata_dict(Dict{Symbol, Any}())
    )
    operators = (
        identity = (; kwargs...) -> qeye(2),
        x = (; kwargs...) -> sigmax(),
        y = (; kwargs...) -> sigmay(),
        z = (; kwargs...) -> sigmaz(),
        p = (; kwargs...) -> sigmap(),
        m = (; kwargs...) -> sigmam(),
    )
    hamiltonian = _component_hamiltonian(hamiltonian, extra_terms)
    return QubitSpec(frequency, dimension, operators, hamiltonian)
end
