Base.show(io::IO, i::Interaction) = print(
    io,
    "Interaction(",
    i.name,
    ", expression=",
    i.expression,
    ", parameters=",
    i.parameters,
    ")",
)

function Base.show(io::IO, ::MIME"text/plain", interaction::Interaction)
    return summary_lines(
        io,
        interaction,
        [
            "Interaction $(interaction.name)",
            "  expression: $(interaction.expression)",
            "  idle parameters: " *
            parameter_values(interaction.parameters),
        ],
    )
end

function mathbody(interaction::Interaction, labels)
    return "H_{" *
           identifier(interaction.name) *
           "} = " *
           render_expression(interaction.expression, true, labels)
end
latex(interaction::Interaction; labels = Dict()) = latex_object(interaction, labels)
Base.show(io::IO, ::MIME"text/latex", interaction::Interaction) = print(io, latex(interaction))
