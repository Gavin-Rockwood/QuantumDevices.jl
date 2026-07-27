"""
    FluxTunableTransmonSpec(EC, EJ1, EJ2; flux = 0, ng = 0, hamiltonian = ..., extra_terms = nothing)

Charge-basis SQUID transmon spec with flux-dependent effective Josephson
energy.
"""
struct FluxTunableTransmonSpec <: AbstractComponentSpec
    EC::DeviceParameter
    EJ1::DeviceParameter
    EJ2::DeviceParameter
    flux::DeviceParameter
    ng::DeviceParameter
    dimension::Dimension
    operators
    hamiltonian::OperatorExpr
end

function FluxTunableTransmonSpec(
        EC::Real,
        EJ1::Real,
        EJ2::Real;
        flux::Real = 0,
        ng::Real = 0,
        metadata = Dict{Symbol,Any}(),
        hamiltonian = missing,
        extra_terms = nothing)
    EC > 0 || error("FluxTunableTransmonSpec EC must be positive.")
    EJ1 > 0 || error("FluxTunableTransmonSpec EJ1 must be positive.")
    EJ2 > 0 || error("FluxTunableTransmonSpec EJ2 must be positive.")

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
    EJ1_parameter = DeviceParameter(
            :EJ1;
            domain = positivedomain(),
            default = Float64(EJ1),
            description = "First junction Josephson energy in GHz.",
            metadata = parameter_metadata,
        )
    EJ2_parameter = DeviceParameter(
            :EJ2;
            domain = positivedomain(),
            default = Float64(EJ2),
            description = "Second junction Josephson energy in GHz.",
            metadata = parameter_metadata,
        )
    flux_parameter = DeviceParameter(
            :flux;
            domain = realdomain(),
            fixed = false,
            default = Float64(flux),
            description = "External flux in flux-quanta units.",
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
        (
            ((param(EJ1_parameter) + param(EJ2_parameter)) * cos(pi * param(flux_parameter)))^2 +
            ((param(EJ1_parameter) - param(EJ2_parameter)) * sin(pi * param(flux_parameter)))^2
        )^(1 / 2) * op(:tunneling) / 2
    resolved_hamiltonian = hamiltonian === missing ? default_hamiltonian : hamiltonian
    return FluxTunableTransmonSpec(
        EC_parameter,
        EJ1_parameter,
        EJ2_parameter,
        flux_parameter,
        ng_parameter,
        dimension,
        operators,
        _component_hamiltonian(resolved_hamiltonian, extra_terms),
    )
end

function FluxTunableTransmonSpec(
        EC::Real,
        EJ::Real;
        d = 0,
        kwargs...
        )
        EJ1 = 0.5*(EJ + d*EJ)
        EJ2 = 0.5*(EJ - d*EJ)
        FluxTunableTransmonSpec(EC, EJ1, EJ2; kwargs...)
end
