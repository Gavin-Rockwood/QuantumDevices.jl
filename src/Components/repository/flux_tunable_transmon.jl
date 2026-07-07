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
    dimension::DeviceParameter
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
        hamiltonian = 4 * param(:EC) * (op(:n) - param(:ng) * op(:identity))^2 -
                      (
                          ((param(:EJ1) + param(:EJ2)) * cos(pi * param(:flux)))^2 +
                          ((param(:EJ1) - param(:EJ2)) * sin(pi * param(:flux)))^2
                      )^(1 / 2) * op(:tunneling) / 2,
        extra_terms = nothing)
    EC > 0 || error("FluxTunableTransmonSpec EC must be positive.")
    EJ1 > 0 || error("FluxTunableTransmonSpec EJ1 must be positive.")
    EJ2 > 0 || error("FluxTunableTransmonSpec EJ2 must be positive.")

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
    return FluxTunableTransmonSpec(
        DeviceParameter(
            :EC;
            domain = positivedomain(),
            default = Float64(EC),
            description = "Charging energy in GHz.",
            metadata = parameter_metadata,
        ),
        DeviceParameter(
            :EJ1;
            domain = positivedomain(),
            default = Float64(EJ1),
            description = "First junction Josephson energy in GHz.",
            metadata = parameter_metadata,
        ),
        DeviceParameter(
            :EJ2;
            domain = positivedomain(),
            default = Float64(EJ2),
            description = "Second junction Josephson energy in GHz.",
            metadata = parameter_metadata,
        ),
        DeviceParameter(
            :flux;
            domain = realdomain(),
            fixed = false,
            default = Float64(flux),
            description = "External flux in flux-quanta units.",
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
