# Recipe summaries and planned-build displays share the graph renderer.
function Base.show(io::IO, recipe::Leaf)
    return print(
        io,
        "Leaf(",
        recipe.component,
        ", dimension=",
        recipe.dimension,
        ", parameters=",
        recipe.overrides,
        ", generator overrides=",
        length(recipe.generators),
        ", control maps=", collect(keys(recipe.control_maps)),
        ")",
    )
end

function Base.show(io::IO, recipe::ResolvedRecipe)
    print(io, "ResolvedRecipe(")
    show(io, recipe.definition)
    return print(
        io,
        ", resolved parameters=",
        recipe.resolved_parameters,
        ", inputs=",
        first.(recipe.inputs),
        ")",
    )
end

function Base.show(io::IO, recipe::Combine)
    return print(
        io,
        "Combine(inputs=",
        recipe.inputs,
        ", interactions=",
        recipe.interactions,
        ", parameters=",
        recipe.overrides,
        ", control maps=", collect(keys(recipe.control_maps)),
        ")",
    )
end

Base.show(io::IO, recipe::CompositeRecipe) = print(
    io,
    "Recipe(components=",
    recipe.components,
    ", inputs=",
    keys(recipe.inputs),
    ", states=",
    Dict(alias => state_description(states) for (alias, states) in recipe.states),
    ", local interactions=",
    [i.name for i in recipe.interactions],
    ", max_dimension=",
    recipe.max_dimension,
    ", control maps=", collect(keys(recipe.control_maps)),
    ")",
)

# Coverage summaries use declarations without resolving or generating matrices.
function coverage_description(recipe)
    recipe isa Leaf && return string(recipe.component)
    recipe isa CompositeRecipe && return join(recipe.components, ", ")
    return "unknown (named upstream inputs)"
end

# Axis coverage stays compact even for sparse retained labels.
function label_axis_description(levels)
    contiguous = last(levels) - first(levels) + 1 == length(levels)
    if !contiguous && length(levels) > 8
        return "$(first(levels)):$(last(levels)) ($(length(levels)) selected levels; gaps)"
    end
    return state_description(levels)
end

# Distinguish Cartesian label products from sparse axis coverage.
function state_description(states)
    isnothing(states) && return "all"
    isempty(states) && return "[]"
    if first(states) isa Tuple
        axes = [sort!(unique(label[index] for label in states)) for index in eachindex(first(states))]
        complete = prod(big(length(axis)) for axis in axes) == length(states)
        ranges = join((string(index, "=", label_axis_description(axis)) for (index, axis) in enumerate(axes)), "; ")
        return "$(length(states)) occupation labels · " *
            (complete ? "Cartesian product" : "sparse selection; axis coverage only") * " · axes " * ranges
    end
    length(states) > 1 &&
        all(i -> states[i] == states[1] + i - 1, eachindex(states)) &&
        return string(first(states), ":", last(states))
    return "[" * join(states, ", ") * "]"
end

# Default recipe inspection delegates to the graph traversal and renderer.
function Base.show(io::IO, ::MIME"text/plain", recipe::Recipe)
    if !(recipe isa Leaf)
        get(io, :compact, false) && return show(io, recipe)
        print(io, recipe isa Combine ? "Combine (compatibility)" : "Recipe",
              " · components: ", coverage_description(recipe), "\n",
              "  control maps: ", list_summary(io,keys(recipe.control_maps)), "\n",
              "Planned build (branches are upstream inputs)\n\n")
        view = graph(recipe)
        lines = String[]
        width = displaysize(io)[2]
        graph_node_rows!(lines, view.root, "", width; root = true)
        graph_lines!(lines, view.root, "", width)
        return print(io, join(lines, "\n"))
    end
    lines = [
        (recipe isa Leaf ? "Leaf" : recipe isa Combine ? "Combine (compatibility)" : "Recipe") *
        " · components: " *
        coverage_description(recipe),
    ]
    if recipe isa Leaf
        push!(
            lines,
            "  initial dimension: " *
            (isnothing(recipe.dimension) ? "default" : string(recipe.dimension)),
        )
        push!(
            lines,
            "  generator overrides: " * list_summary(io, [r.name for r in keys(recipe.generators)]),
        )
    end
    push!(lines, "  control maps: " * list_summary(io,keys(recipe.control_maps)))
    push!(lines, "  declared parameter overrides: " * parameter_values(recipe.overrides))
    return summary_lines(io, recipe, lines)
end

# Resolved provenance summarizes captured values separately from declared overrides.
function Base.show(io::IO, ::MIME"text/plain", recipe::ResolvedRecipe)
    return summary_lines(
        io,
        recipe,
        [
            "ResolvedRecipe · components: " * coverage_description(recipe.definition),
            "  resolved build values: " * parameter_values(recipe.resolved_parameters),
            "  immediate input aliases: " * list_summary(io, first.(recipe.inputs); ordered = true),
            "  inspect .definition for declared overrides or graph(recipe) for upstream inputs",
        ],
    )
end


Base.show(io::IO, binding::InputBinding) =
    print(io, "InputBinding(", binding.name, ", owners=", binding.owners, ")")
Base.show(io::IO, operators::InputOperators) =
    print(io, "InputOperators(", getfield(operators, :binding).name, ", names=", propertynames(operators), ")")

