"""
    ResonatorSpec(frequency; dimension = 10, hamiltonian = param(:frequency) * op(:n))

Truncated harmonic resonator with ladder, number, and quadrature operators.
"""
struct ResonatorSpec <: AbstractComponentSpec
    frequency::DeviceParameter
    dimension::DeviceParameter
    operators
    hamiltonian::OperatorExpr
end

function ResonatorSpec(
        frequency::Real;
        dimension::Integer = 10,
        metadata = Dict{Symbol,Any}(),
        hamiltonian = param(:frequency) * op(:n),
        extra_terms = nothing)
    frequency > 0 || error("ResonatorSpec frequency must be positive.")
    dimension > 0 || error("ResonatorSpec dimension must be positive.")

    parameter_metadata = _metadata_dict(metadata)
    operators = (
        a = (; dimension, kwargs...) -> destroy(dimension),
        adag = (; dimension, kwargs...) -> create(dimension),
        n = (; dimension, kwargs...) -> num(dimension),
        number = (; dimension, kwargs...) -> num(dimension),
        q = (; dimension, kwargs...) -> destroy(dimension) + create(dimension),
        p = (; dimension, kwargs...) -> -1im * (destroy(dimension) - create(dimension)),
    )
    return ResonatorSpec(
        DeviceParameter(
            :frequency;
            domain = positivedomain(),
            default = Float64(frequency),
            description = "Resonator frequency in GHz.",
            metadata = parameter_metadata,
        ),
        DeviceParameter(
            :dimension;
            domain = integerrange(1, typemax(Int)),
            fixed = true,
            default = Int(dimension),
            description = "Resonator Fock-space dimension.",
            metadata = parameter_metadata,
        ),
        operators,
        _component_hamiltonian(hamiltonian, extra_terms),
    )
end
