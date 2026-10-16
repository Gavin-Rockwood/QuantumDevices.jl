"""
    Component(name, parameters, operators, hamiltonian, type, dimension)
    Component(component; operators=component.operators, hamiltonian=component.hamiltonian)

A component with a string name, named scalar parameters, named local quantum
operators, a symbolic Hamiltonian, a storage type tag, and a parent dimension.
All operators must be square matrices of that dimension when used in a model.
Use the built-in `make_*` constructors for basis consistency; this constructor
also supports custom components. The copy form preserves physical basis metadata
when extending operators in the same coordinate basis. Component names must be unique within a model.
"""
struct Component
    name::AbstractString
    parameters::NamedTuple
    operators::NamedTuple
    hamiltonian
    type::String
    dimension
    basis
end

# Each built-in component registers its keyword constructor. Stored parameters
# can then be splatted into it without depending on their serialized order.
const COMPONENT_CONSTRUCTORS = Dict{String,Function}()

# The map sends component coordinates into a declared common physical basis.
struct _ComponentBasis
    id
    vectors
end
Component(name,parameters,operators,hamiltonian,type,dimension;basis=nothing) =
    Component(name,parameters,operators,hamiltonian,type,dimension,basis)
Component(component::Component;operators=component.operators,hamiltonian=component.hamiltonian) =
    Component(component.name,component.parameters,operators,hamiltonian,
        component.type,component.dimension;basis=component.basis)
