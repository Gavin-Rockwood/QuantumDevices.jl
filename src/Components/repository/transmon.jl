"""
    TransmonSpec(EC, EJ; ng = 0, hamiltonian = ..., extra_terms = nothing)

Charge-basis transmon spec whose Hamiltonian is built from identity, centered
charge (`:n`), and tunneling operators.
"""
struct TransmonSpec <: AbstractComponentSpec
    EC::DeviceParameter
    EJ::DeviceParameter
    ng::DeviceParameter
    dimension::DeviceParameter
    operators
    hamiltonian::OperatorExpr
end

function TransmonSpec(
        EC::Real,
        EJ::Real;
        ng::Real = 0,
        metadata = Dict{Symbol,Any}(),
        hamiltonian = 4 * param(:EC) * (op(:n) - param(:ng) * op(:identity))^2 -
                      param(:EJ) * op(:tunneling) / 2,
        extra_terms = nothing)
    EC > 0 || error("TransmonSpec EC must be positive.")
    EJ > 0 || error("TransmonSpec EJ must be positive.")
    
    parameter_metadata = _metadata_dict(metadata)
    dimension_parameter = DeviceParameter(
        :dimension;
        domain = anydomain(),
        fixed = true,
        default = Inf,
        description = "Maximum charge-basis dimension.",
        metadata = parameter_metadata,
    )
    operators = (
        identity = (; dimension, kwargs...) ->
            qeye(_transmon_dimension(dimension)),
        n = (; dimension, kwargs...) -> begin
            dim = _transmon_dimension(dimension)
            num(dim) - (dim ÷ 2) * qeye(dim)
        end,
        tunneling = (; dimension, kwargs...) ->
            tunneling(_transmon_dimension(dimension)),
    )
    return TransmonSpec(
        DeviceParameter(
            :EC;
            domain = positivedomain(),
            default = Float64(EC),
            description = "Charging energy in GHz.",
            metadata = parameter_metadata,
        ),
        DeviceParameter(
            :EJ;
            domain = positivedomain(),
            default = Float64(EJ),
            description = "Josephson energy in GHz.",
            metadata = parameter_metadata,
        ),
        DeviceParameter(
            :ng;
            domain = realdomain(),
            default = Float64(ng),
            description = "Offset charge.",
            metadata = parameter_metadata,
        ),
        dimension_parameter,
        operators,
        _component_hamiltonian(hamiltonian, extra_terms),
    )
end

"""
    _transmon_dimension(dimension)

Validate and normalize the finite odd charge-basis dimension required for
Transmon operators.
"""
function _transmon_dimension(dimension)
    isfinite(dimension) ||
        error("TransmonSpec requires an explicit finite numerical dimension.")
    dimension isa Integer && dimension > 0 && isodd(dimension) ||
        error("Transmon numerical dimension must be a positive odd integer.")
    return Int(dimension)
end
