"""
    TransmonSpec(EC, EJ; ng = 0, hamiltonian = ..., extra_terms = nothing)

Charge-basis transmon spec whose Hamiltonian is built from identity, centered
charge (`:n`), and tunneling operators.
"""
struct TransmonSpec <: AbstractComponentSpec
    EC::DeviceParameter
    EJ::DeviceParameter
    ng::DeviceParameter
    dimension::Dimension
    operators
    hamiltonian::OperatorExpr
end

function TransmonSpec(
        EC::Real,
        EJ::Real;
        ng::Real = 0,
        metadata = Dict{Symbol,Any}(),
        hamiltonian = missing,
        extra_terms = nothing)
    EC > 0 || error("TransmonSpec EC must be positive.")
    EJ > 0 || error("TransmonSpec EJ must be positive.")
    
    parameter_metadata = _metadata_dict(metadata)
    dimension = Dimension(Inf)
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
    EC_parameter = DeviceParameter(
            :EC;
            domain = positivedomain(),
            default = Float64(EC),
            description = "Charging energy in GHz.",
            metadata = parameter_metadata,
        )
    EJ_parameter = DeviceParameter(
            :EJ;
            domain = positivedomain(),
            default = Float64(EJ),
            description = "Josephson energy in GHz.",
            metadata = parameter_metadata,
        )
    ng_parameter = DeviceParameter(
            :ng;
            domain = realdomain(),
            default = Float64(ng),
            description = "Offset charge.",
            metadata = parameter_metadata,
        )
    default_hamiltonian =
        4 * param(EC_parameter) * (op(:n) - param(ng_parameter) * op(:identity))^2 -
        param(EJ_parameter) * op(:tunneling) / 2
    resolved_hamiltonian = hamiltonian === missing ? default_hamiltonian : hamiltonian
    return TransmonSpec(
        EC_parameter,
        EJ_parameter,
        ng_parameter,
        dimension,
        operators,
        _component_hamiltonian(resolved_hamiltonian, extra_terms),
    )
end

"""
    _transmon_dimension(dimension)

Validate and normalize the finite odd charge-basis dimension required for
Transmon operators.
"""
function _transmon_dimension(dimension)
    dimension isa Tuple && (dimension = only(dimension))
    isfinite(dimension) ||
        error("TransmonSpec requires an explicit finite numerical dimension.")
    dimension isa Integer && dimension > 0 && isodd(dimension) ||
        error("Transmon numerical dimension must be a positive odd integer.")
    return Int(dimension)
end
