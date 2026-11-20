"""Construct an oscillator with `frequency * n`; `subscript` selects its variable suffix."""
function make_resonator(name; frequency = 3.0, dimension = 4, label = string(name), metadata = (;), subscript = missing)
    variables, parameters = component_parameters(name; subscript, ν=frequency)
    ν = variables["ν"]
    generators = Dict(
        "a" => OperatorGenerator(BuiltinGenerator(:destroy), dimension),
        "n" => OperatorGenerator(BuiltinGenerator(:number), dimension),
    )
    operators, operator_generators = component_operators(name,
        generators)
    return Component(name; variables, parameters, operators, operator_generators,
        hamiltonian=ν * op(name, :n), type="resonator", label, metadata)
end
