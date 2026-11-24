function pauli_dimension(dimension)
    validate_operator_dimension(dimension) == 2 || throw(ArgumentError("Pauli basis needs dimension 2"))
    return 2
end
pauli_x_operator(dimension) = (pauli_dimension(dimension); sigmax())
pauli_y_operator(dimension) = (pauli_dimension(dimension); sigmay())
pauli_z_operator(dimension) = (pauli_dimension(dimension); sigmaz())
merge!(BUILTIN_OPERATOR_FUNCTIONS, Dict{String,Function}("pauli_x" => pauli_x_operator, "pauli_y" => pauli_y_operator, "pauli_z" => pauli_z_operator))

"""
    make_qubit(name::Symbol, frequency; label=missing, metadata=(;), subscript=missing)
    make_qubit(name::Symbol; frequency=5.0, kwargs...)
    make_qubit(; name, frequency=5.0, kwargs...)

Construct a Pauli qubit with `frequency * z / 2`; `subscript` selects its variable suffix.
"""
function make_qubit(name::Symbol, frequency; label = missing, metadata = (;), subscript = missing)
    variables, parameters = component_parameters(name; subscript, ν=frequency)
    ν = variables["ν"]
    operators, operator_generators = component_operators(name,
        Dict{String,Function}("x" => pauli_x_operator, "y" => pauli_y_operator, "z" => pauli_z_operator))
    return Component(name; variables, parameters, operators, operator_generators,
        dimension=2, hamiltonian=0.5 * ν * op(name, :z), type="qubit",
        label=ismissing(label) ? string(name) : label, metadata)
end

make_qubit(name::Symbol; frequency = 5.0, kwargs...) =
    make_qubit(name, frequency; kwargs...)
make_qubit(; name, frequency = 5.0, kwargs...) =
    make_qubit(name, frequency; kwargs...)
