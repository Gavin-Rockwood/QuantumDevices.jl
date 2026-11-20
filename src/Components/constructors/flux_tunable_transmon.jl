"""
    make_flux_tunable_transmon(name::Symbol, EC, EJmax;
        d=0.0, ng=0.0, flux=0.0, dimension=101, label=missing, metadata=(;), subscript=missing)
    make_flux_tunable_transmon(name::Symbol; EC, EJmax, kwargs...)
    make_flux_tunable_transmon(; name, EC, EJmax, kwargs...)

Construct a split-junction transmon with effective Josephson energy
`EJmax * sqrt(cospi(flux)^2 + d^2*sinpi(flux)^2)`.
`flux` is external flux `Φext/Φ0`, and `d` is junction asymmetry.
All five physical inputs remain independent symbolic parameters. Charge-basis
generators, dimensions, and energy units follow [`make_transmon`](@ref).
"""
function make_flux_tunable_transmon(
    name::Symbol,
    EC,
    EJmax;
    d = 0.0,
    ng = 0.0,
    flux = 0.0,
    dimension = 101,
    label = missing,
    metadata = (;),
    subscript = missing,
)
    variables, parameters = component_parameters(name; subscript, EC, EJmax, d, ng, flux)
    effective_EJ =
        variables["EJmax"] * sqrt(
            cospi(variables["flux"]) ^ 2 +
            variables["d"] ^ 2 * sinpi(variables["flux"]) ^ 2,
        )
    hamiltonian =
        4 * variables["EC"] * (op(name, :charge) - variables["ng"]) ^ 2 -
        effective_EJ * op(name, :jump) / 2
    operators, operator_generators = component_operators(name,
        charge_basis_generators(name, dimension))
    return Component(name; variables, parameters, operators, operator_generators, hamiltonian,
        type="flux-tunable transmon", label=ismissing(label) ? string(name) : label, metadata)
end

make_flux_tunable_transmon(name::Symbol; EC, EJmax, kwargs...) =
    make_flux_tunable_transmon(name, EC, EJmax; kwargs...)
make_flux_tunable_transmon(; name, EC, EJmax, kwargs...) =
    make_flux_tunable_transmon(name, EC, EJmax; kwargs...)
