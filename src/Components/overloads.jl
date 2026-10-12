"""
    numerical(component::Component, expression::Sym; scalar=1, dense=false)

Evaluate a symbolic expression using the component's local operators and
parameters in its parent Hilbert space. Use local names such as `op(:charge)`.
For example, `numerical(component, component.hamiltonian)` evaluates its Hamiltonian.
Unresolved parameters produce a `QobjEvo`, as with the other symbolic overloads.
`dense=true` materializes all operator matrices as dense matrices at construction.
`scalar` is folded into operator matrices and constant coefficients at construction.
"""
numerical(component::Component, expression::Sym; scalar::Number=1, dense::Bool=false) =
    numerical(expression, component.operators, component.parameters; scalar, dense)

function _replace_property(component::Component, key::Symbol, value)
    return _rebuild_component(component, Dict(key => value))
end

function _rebuild_component(component::Component, changes)
    all(key -> key in (:name, :parameters, :dimension), keys(changes)) ||
        throw(ArgumentError("Component inputs are derived or unknown"))
    name = get(changes, :name, component.name)
    p = get(changes, :parameters, component.parameters)
    dimension = get(changes, :dimension, component.dimension)
    name === component.name && p === component.parameters && dimension === component.dimension && return component
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
