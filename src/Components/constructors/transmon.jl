# Shared charge-basis physics for fixed and flux-tunable transmons.
function charge_basis_dimension(dimension)
    dimension = validate_operator_dimension(dimension)
    isodd(dimension) || throw(ArgumentError("charge-basis dimension must be a positive odd state count (2*ncut + 1)"))
    return dimension
end

charge_operator(dimension) = num(charge_basis_dimension(dimension)) - (dimension-1) ÷ 2
jump_operator(dimension) = tunneling(charge_basis_dimension(dimension), 1)
merge!(BUILTIN_OPERATOR_FUNCTIONS, Dict{String,Function}("charge" => charge_operator, "jump" => jump_operator))

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
    dimension = charge_basis_dimension(dimension)
    variables, parameters = component_parameters(name; subscript, EC, EJ, ng)
    hamiltonian =
        4 * variables["EC"] * (op(name, :charge) - variables["ng"]) ^ 2 -
        variables["EJ"] * op(name, :jump) / 2
    operators, operator_generators = component_operators(name,
        Dict{String,Function}("charge" => charge_operator, "jump" => jump_operator))
    return Component(name; variables, parameters, operators, operator_generators, hamiltonian, dimension,
        type="transmon", label=ismissing(label) ? string(name) : label, metadata)
end

make_transmon(name::Symbol; EC, EJ, kwargs...) = make_transmon(name, EC, EJ; kwargs...)
make_transmon(; name, EC, EJ, kwargs...) = make_transmon(name, EC, EJ; kwargs...)
