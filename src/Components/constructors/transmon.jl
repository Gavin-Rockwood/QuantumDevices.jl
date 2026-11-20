function charge_basis_dimension(dimension)
    dimension isa Integer && dimension > 0 && isodd(dimension) || throw(
        ArgumentError("charge-basis dimension must be a positive odd state count (2*ncut + 1)"),
    )
    return dimension
end

function charge_basis_generators(name::Symbol, dimension)
    charge_basis_dimension(dimension)
    return Dict(
        "charge" => OperatorGenerator(
            BuiltinGenerator(:charge),
            dimension,
        ),
        "jump" =>
            OperatorGenerator(BuiltinGenerator(:jump), dimension),
    )
end

"""
    make_transmon(name::Symbol, EC, EJ; ng=0.0, dimension=101, label=missing, metadata=(;), subscript=missing)
    make_transmon(name::Symbol; EC, EJ, kwargs...)
    make_transmon(; name, EC, EJ, kwargs...)

Construct the symbolic charge-basis Hamiltonian `4EC*(charge-ng)^2 - EJ*jump/2`.
`EC` and `EJ` use frequency units `E/h`; `ng` is dimensionless offset charge.
The `charge` and `jump` operator generators return QuantumToolbox operators,
with charge states `-ncut:ncut` and nearest-neighbor tunneling respectively.
`dimension = 2*ncut + 1` is a positive odd state count, defaulting to 101.
No operators are generated or diagonalized until a numerical model is built.
Increase the charge dimension to check convergence of retained energy levels.
"""
function make_transmon(
    name::Symbol,
    EC,
    EJ;
    ng = 0.0,
    dimension = 101,
    label = missing,
    metadata = (;),
    subscript = missing,
)
    variables, parameters = component_parameters(name; subscript, EC, EJ, ng)
    hamiltonian =
        4 * variables["EC"] * (op(name, :charge) - variables["ng"]) ^ 2 -
        variables["EJ"] * op(name, :jump) / 2
    operators, operator_generators = component_operators(name,
        charge_basis_generators(name, dimension))
    return Component(name; variables, parameters, operators, operator_generators, hamiltonian,
        type="transmon", label=ismissing(label) ? string(name) : label, metadata)
end

make_transmon(name::Symbol; EC, EJ, kwargs...) = make_transmon(name, EC, EJ; kwargs...)
make_transmon(; name, EC, EJ, kwargs...) = make_transmon(name, EC, EJ; kwargs...)
