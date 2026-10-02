"""
    make_resonator(name, frequency, dimension)
    make_resonator(; name, frequency, dimension)

Construct a harmonic resonator in its Fock basis with
`H = frequency * a†a`. `frequency` is in cycles per unit time
(e.g. GHz for time in ns); evolution uses `2π * H`. The positive integer `dimension` is the
number of Fock states retained in the component's parent space.

The local operators are `a` (annihilation), `adag` (creation), and `n` (number).
The vacuum energy is omitted; it would only add a constant to the Hamiltonian.
Check convergence by increasing `dimension` when higher occupations matter.
"""
function make_resonator(name, frequency, dimension)
    dimension isa Integer && dimension > 0 ||
        throw(ArgumentError("Resonator dimension must be a positive integer"))

    a = destroy(dimension)
    operators = (a=a, adag=a', n=num(dimension))
    parameters = (frequency=frequency,)
    hamiltonian = param(:frequency) * op(:n)
    return Component(name, parameters, operators, hamiltonian, "Resonator", dimension)
end

make_resonator(; name, frequency, dimension) = make_resonator(name, frequency, dimension)

COMPONENT_CONSTRUCTORS["Resonator"] = make_resonator
