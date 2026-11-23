"""
    make_qubit(name::Symbol, frequency; label=missing, metadata=(;), subscript=missing)
    make_qubit(name::Symbol; frequency=5.0, kwargs...)
    make_qubit(; name, frequency=5.0, kwargs...)

Construct a Pauli qubit with `frequency * z / 2`; `subscript` selects its variable suffix.
"""
function make_qubit(name::Symbol, frequency; label = missing, metadata = (;), subscript = missing)
    variables, parameters = component_parameters(name; subscript, ν=frequency)
    ν = variables["ν"]
    paulis = Dict(:x => sigmax, :y => sigmay, :z => sigmaz)
    generators = Dict(
        string(local_name) => OperatorGenerator(
            BuiltinGenerator(Symbol("pauli_",local_name)),
            2,
        ) for (local_name, matrix) in paulis
    )
    operators, operator_generators = component_operators(name,
        generators)
    return Component(name; variables, parameters, operators, operator_generators,
        hamiltonian=0.5 * ν * op(name, :z), type="qubit",
        label=ismissing(label) ? string(name) : label, metadata)
end

make_qubit(name::Symbol; frequency = 5.0, kwargs...) =
    make_qubit(name, frequency; kwargs...)
make_qubit(; name, frequency = 5.0, kwargs...) =
    make_qubit(name, frequency; kwargs...)
