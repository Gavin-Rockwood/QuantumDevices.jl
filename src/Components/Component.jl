"""
    Component(name, parameters, operators, hamiltonian, type, dimension)

A component with a string name, named scalar parameters, named local quantum
operators, a symbolic Hamiltonian, a storage type tag, and a parent dimension.
All operators must be square matrices of that dimension when used in a model.
Use the built-in `make_*` constructors for basis consistency; this constructor
also supports custom components. Component names must be unique within a model.
"""
struct Component
    name::AbstractString
    parameters::NamedTuple
    operators::NamedTuple
    hamiltonian
    type::String
    dimension
end

# Each built-in component registers its keyword constructor. Stored parameters
# can then be splatted into it without depending on their serialized order.
const COMPONENT_CONSTRUCTORS = Dict{String,Function}()
