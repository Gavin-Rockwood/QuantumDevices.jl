"""
    ResonatorSpec(frequency; dimension = 10, hamiltonian = missing)

Truncated harmonic resonator with ladder, number, and quadrature operators.
"""
struct ResonatorSpec <: AbstractComponentSpec
    frequency::DeviceParameter
    dimension::Dimension
    operators
    hamiltonian::OperatorExpr
end

function ResonatorSpec(
        frequency::Real;
        dimension::Integer = 10,
        metadata = Dict{Symbol,Any}(),
        hamiltonian = missing,
        extra_terms = nothing)
    frequency > 0 || error("ResonatorSpec frequency must be positive.")
    dimension > 0 || error("ResonatorSpec dimension must be positive.")

    parameter_metadata = _metadata_dict(metadata)
    operators = (
        a = (; dimension, kwargs...) -> destroy(only(dimension)),
        adag = (; dimension, kwargs...) -> create(only(dimension)),
        n = (; dimension, kwargs...) -> num(only(dimension)),
        number = (; dimension, kwargs...) -> num(only(dimension)),
        q = (; dimension, kwargs...) -> destroy(only(dimension)) + create(only(dimension)),
        p = (; dimension, kwargs...) -> -1im * (destroy(only(dimension)) - create(only(dimension))),
    )
    frequency_parameter = DeviceParameter(
            :frequency;
            domain = positivedomain(),
            default = Float64(frequency),
            description = "Resonator frequency in GHz.",
            metadata = parameter_metadata,
        )
    dimension_config = Dimension(dimension)
    default_hamiltonian = param(frequency_parameter) * op(:n)
    resolved_hamiltonian = hamiltonian === missing ? default_hamiltonian : hamiltonian
    return ResonatorSpec(
        frequency_parameter,
        dimension_config,
        operators,
        _component_hamiltonian(resolved_hamiltonian, extra_terms),
    )
end
