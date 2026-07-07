"""
    GenericSpec(spectrum; hamiltonian = op(:energy), extra_terms = nothing)

Finite-dimensional component with a diagonal energy spectrum and identity
operator.
"""
struct GenericSpec <: AbstractComponentSpec
    spectrum::Vector{DeviceParameter}
    dimension::DeviceParameter
    operators
    hamiltonian::OperatorExpr
end

function GenericSpec(vals;
    domain = realdomain(),
    description = "System Spectrum",
    metadata = Dict{Symbol, Any}(),
    hamiltonian = op(:energy),
    extra_terms = nothing
    )
    spectrum = [DeviceParameter(
        ParamPath(:spectrum),
        domain,
        true,
        true,
        val,
        description,
        _metadata_dict(metadata)
    ) for val in vals]

    dimension = DeviceParameter(
        ParamPath(:dimension),
        integerrange(length(vals),length(vals)),
        true,
        true,
        length(vals),
        "System Dimension",
        _metadata_dict(Dict{Symbol, Any}())
    )
    spectrum_values = [parameter.default for parameter in spectrum]
    operators = (
        identity = (; dimension, kwargs...) -> qeye(dimension),
        energy = (; dimension, kwargs...) -> begin
            dimension == length(spectrum_values) ||
                error("GenericSpec dimension $dimension does not match spectrum length $(length(spectrum_values)).")
            QuantumObject(Diagonal(ComplexF64.(spectrum_values)))
        end,
    )
    return GenericSpec(
        spectrum,
        dimension,
        operators,
        _component_hamiltonian(hamiltonian, extra_terms),
    )
end
