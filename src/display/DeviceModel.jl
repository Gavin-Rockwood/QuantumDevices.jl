function Base.show(io::IO, model::DeviceModel)
    print(io, "DeviceModel(", length(model.components), " components, dimension=",
        size(model.H, 1), ")")
end

function Base.show(io::IO, ::MIME"text/plain", model::DeviceModel)
    show(io, model)
    limited = IOContext(io, :compact => true, :limit => true)
    print(io, "\n  Components: ")
    show(limited, [c.name for c in model.components])
    print(io, "\n  Parameters: ")
    show(limited, model.parameters)
    print(io, "\n  Operators: ", length(model.operators))
    print(io, "\n  Compiled terms: ", length(model.compiled_terms))
    print(io, "\n  Gates: ")
    show(limited, collect(keys(model.gates)))
end
