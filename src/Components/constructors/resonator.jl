annihilation_operator(dimension) = destroy(dimension)
number_operator(dimension) = num(dimension)
merge!(BUILTIN_OPERATOR_FUNCTIONS, Dict{String,Function}("destroy" => annihilation_operator, "number" => number_operator))

"""
    make_resonator(name, frequency; kwargs...)
    make_resonator(name; frequency=3.0, kwargs...)
    make_resonator(; name, frequency=3.0, kwargs...)

Construct an oscillator with `frequency * n`; `subscript` selects its variable suffix.
"""
function make_resonator(name::Symbol, frequency; dimension = 4, label = missing, metadata = (;), subscript = missing)
    label = ismissing(label) ? string(name) : label
    variables, parameters = component_parameters(name; subscript, ν=frequency)
    ν = variables["ν"]
    operators, operator_generators = component_operators(name,
        Dict{String,Function}("a" => annihilation_operator, "n" => number_operator))
    return Component(name; variables, parameters, operators, operator_generators,
        dimension, hamiltonian=ν * op(name, :n), type="resonator", label, metadata)
end

make_resonator(name::Symbol; frequency = 3.0, kwargs...) =
    make_resonator(name, frequency; kwargs...)

make_resonator(; name, frequency = 3.0, kwargs...) =
    make_resonator(name, frequency; kwargs...)
