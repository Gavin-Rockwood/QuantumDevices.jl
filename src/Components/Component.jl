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

function _replace_property(component::Component, key::Symbol, value)
    key in (:name, :parameters, :dimension) ||
        throw(ArgumentError("Component.$key is derived or unknown"))
    name = key === :name ? value : component.name
    p = key === :parameters ? value : component.parameters
    dimension = key === :dimension ? value : component.dimension
    constructor = get(COMPONENT_CONSTRUCTORS, component.type, nothing)
    constructor === nothing &&
        return Component(name, p, component.operators, component.hamiltonian, component.type, dimension)
    rebuilt = constructor(; name, dimension, p...)
    for parameter in keys(p)
        if p[parameter] != component.parameters[parameter] &&
                p[parameter] != rebuilt.parameters[parameter]
            throw(ArgumentError("$parameter is derived from physical parameters"))
        end
    end
    return rebuilt
end
