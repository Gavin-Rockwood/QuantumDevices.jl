function Base.getindex(generators::Dict{OperatorRef, OperatorGenerator}, reference::Sym)
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
    operators = sort!(collect(component.operator_generators); by = p -> first(p))
    push!(
        lines,
        "  operators (default dimensions): " * list_summary(
            io,
            ["$(ref) ($(g.default_dimension))" for (ref, g) in operators];
            ordered = true,
        ),
    )
    push!(lines, "  inspect .hamiltonian for the symbolic Hamiltonian")
    return summary_lines(io, component, lines)
end

function (g::OperatorGenerator)(dimension = g.default_dimension; kwargs...)
    dimension isa Integer && dimension > 0 ||
        throw(ArgumentError("dimension must be a positive state count"))
    return g.f(dimension; kwargs...)
end
