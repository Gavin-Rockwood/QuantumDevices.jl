"""
    numerical(component::Component, expression::Sym)

Evaluate a symbolic expression using the component's local operators and
parameters in its parent Hilbert space. Use local names such as `op(:charge)`.
For example, `numerical(component, component.hamiltonian)` evaluates its Hamiltonian.
Unresolved parameters produce a `QobjEvo`, as with the other symbolic overloads.
"""
numerical(component::Component, expression::Sym) =
    numerical(expression, component.operators, component.parameters)

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

function Base.show(io::IO, component::Component)
    print(io, "Component(")
    show(io, component.name)
    print(io, ", ", component.type, ", dimension=", component.dimension, ")")
end

function Base.show(io::IO, ::MIME"text/plain", component::Component)
    show(io, component)
    limited = IOContext(io, :compact => true, :limit => true)
    print(io, "\n  Parameters: ")
    show(limited, component.parameters)
    print(io, "\n  Operators: ")
    show(limited, keys(component.operators))
    print(io, "\n  Hamiltonian: ")
    show(limited, component.hamiltonian)
end
