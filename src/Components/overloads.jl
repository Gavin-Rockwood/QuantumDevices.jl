function Base.getindex(generators::Dict{OperatorRef, Function}, reference::Sym)
    reference.expr[1] == :op || throw(ArgumentError("Expected operator reference"))
    return generators[reference.expr[2]]
end

Base.show(io::IO, c::Component) = print(
    io,
    "Component(",
    c.name,
    ", type=",
    repr(c.type),
    ", parameters=",
    c.parameters,
    ", operators=",
    length(c.operators),
    ")",
)

function Base.show(io::IO, ::MIME"text/plain", component::Component)
    lines = ["Component $(component.name) · $(component.type)"]
    !isnothing(component.label) &&
        component.label != string(component.name) &&
        push!(lines, "  label: $(component.label)")
    push!(lines, "  variables: " * join((string(alias," => ",variable) for (alias,variable) in pairs(component.variables)), ", "))
    push!(
        lines,
        "  idle parameters: " *
        parameter_values(component_parameter_values(component)),
    )
    push!(lines, "  default dimension: " * (isnothing(component.dimension) ? "unspecified" : string(component.dimension)))
    operators = sort!(collect(component.operator_generators); by = p -> first(p))
    push!(
        lines,
        "  operators: " * list_summary(
            io,
            [string(ref) for (ref, _) in operators];
            ordered = true,
        ),
    )
    push!(lines, "  inspect .hamiltonian for the symbolic Hamiltonian")
    return summary_lines(io, component, lines)
end
