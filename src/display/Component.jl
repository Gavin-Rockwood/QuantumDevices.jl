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
