"""A read-only description of one upstream build occurrence, independent of layout."""
struct RecipeGraphNode
    path::Tuple
    kind::Symbol # leaf, combine, or terminal unresolved/cyclic input
    name::Union{Nothing, Symbol}
    component::Union{Nothing, Symbol}
    message::String
    dimension::Union{Nothing, BigInt}
    details::Vector{String}
    edges::Vector{NamedTuple} # ordered (; name, states, child), independent of display text
end

"""A displayable recipe construction tree. Create with [`graph`](@ref)."""
struct RecipeGraph
    root::RecipeGraphNode
end

# Infer dimensions from declarations only; inspection never invokes a generator.
function inspection_dimension(device, leaf::Leaf)
    !isnothing(leaf.dimension) && return big(leaf.dimension)
    isnothing(device) && return nothing
    haskey(device.components, leaf.component) || return nothing
    dimension = device.components[leaf.component].dimension
    return isnothing(dimension) ? nothing : big(dimension)
end

# Describe planned execution using the same reduction classification as model assembly.
function planned_stage_details(definition, sources)
    operation = recipe_stage_operation(definition, sources)
    operation === :initialize && return ["Plan: diagonalize component; initialize occupation labels (no state tracking)"]
    if operation === :projection
        return ["Plan: pure reduction of upstream dressed states",
                "Eigensystem: reuse selected eigenstates; no new eigensolve",
                definition.dressing ?
                "Tracking: inherit upstream labels and diagnostics; no new state tracking" :
                "Tracking: disabled; discard upstream labels and diagnostics"]
    elseif !definition.dressing
        plan = isempty(definition.interactions) ?
            "Plan: assemble diagonal product Hamiltonian; no new eigensolve" :
            "Plan: assemble Hamiltonian; diagonalize for nonzero new interaction"
        return [plan, "Tracking: disabled; no occupation labels"]
    elseif operation === :unresolved
        return ["Plan: resolve named upstream to determine projection or dressing",
                "Tracking: depends on resolved upstream; supply graph(device, recipe)"]
    end
    options = definition.dressing_options
    samples = isempty(definition.interactions) ? "2 samples (zero new interaction)" :
        "up to $(options.nsteps) samples"
    plan = isempty(definition.interactions) ?
        "Plan: assemble diagonal product Hamiltonian; no new eigensolve" :
        "Plan: assemble Hamiltonian; diagonalize for nonzero new interaction"
    return [plan,
            "Tracking: perform dressed-state tracking ($samples, :$(options.method), thread_count=$(options.thread_count))"]
end
# Unresolved and cyclic branches remain explanatory terminal nodes.
unknown_graph(path, title) =
    RecipeGraphNode(path, :terminal, nothing, nothing, title, nothing, String[], NamedTuple[])
# Expand each input occurrence, stopping only active-path cycles.
function describe_graph(device, source, path = (); active = Any[])
    savedname = source isa Symbol ? source : nothing
    if source isa Symbol
        isnothing(device) &&
            return unknown_graph(path, "unresolved :$source (supply graph(device, recipe))")
        haskey(device.recipes, source) ||
            return unknown_graph(path, "unresolved :$source (missing registered recipe)")
        source = device.recipes[source]
    end
    any(x -> x === source, active) && return unknown_graph(
        path,
        "cycle detected" * (isnothing(savedname) ? "" : " at :$savedname"),
    )
    resolved = source isa ResolvedRecipe
    definition = resolved ? source.definition : source
    definition isa Recipe || return unknown_graph(
        path,
        "unsupported input $(typeof(source)); expected Recipe or registered name",
    )
    nextactive = vcat(active, Any[source])
    details = String[]
    parameters = resolved ? source.resolved_parameters : definition.overrides
    !isempty(parameters) && push!(
        details,
        (resolved ? "resolved build values: " : "declared overrides: ") *
        parameter_values(parameters),
    )
    edges = NamedTuple[]
    if definition isa Leaf
        dimension = inspection_dimension(device, definition)
        append!(details, planned_stage_details(definition, []))
        return RecipeGraphNode(path, :leaf, savedname, definition.component, "", dimension, details, edges)
    elseif !(definition isa Union{CompositeRecipe, Combine})
        return unknown_graph(path, "unsupported recipe $(typeof(definition))")
    end
    inputs = recipe_inputs(source)
    aliases = [input.name for input in inputs]
    sources = [input.recipe for input in inputs]
    selections = [input.states for input in inputs]
    if !isempty(definition.interactions)
        heading = definition isa CompositeRecipe ? "local interactions: " : "registered physical interactions: "
        push!(details, heading * join(interaction_name.(definition.interactions), ", "))
    end
    # Use resolved definitions for operation inspection without doing numerical work.
    operation_sources = [upstream isa Symbol && !isnothing(device) && haskey(device.recipes, upstream) ?
                         device.recipes[upstream] : upstream for upstream in sources]
    dimensions = Union{Nothing, BigInt}[]
    for (index, (alias, upstream, selection)) in enumerate(zip(aliases, sources, selections))
        child = describe_graph(device, upstream, (path..., index); active = nextactive)
        push!(edges, (; name = alias, states = selection, child))
        push!(dimensions, isnothing(selection) ? child.dimension : big(length(selection)))
    end
    dimension = isempty(dimensions) || any(isnothing, dimensions) ? nothing : prod(dimensions)
    append!(details, planned_stage_details(definition, operation_sources))
    push!(
        details,
        "max_dimension: " *
        (definition.max_dimension == 0 ? "unlimited" : string(definition.max_dimension)),
    )
    return RecipeGraphNode(path, :combine, savedname, nothing, "", dimension, details, edges)
end

"""
    graph(recipe)
    graph(device, recipe_or_name)

Inspect the entire upstream construction without generating operators or enforcing
build limits. Embedded objects expand directly; supply a device to resolve saved
names and component default dimensions. ResolvedRecipe snapshots need no device.
Shared inputs expand per occurrence; only active-path cycles stop a branch.
"""
graph(source::Union{Recipe, ResolvedRecipe}) = RecipeGraph(describe_graph(nothing, source))
graph(device::QuantumDevice, source::Union{Recipe, ResolvedRecipe, Symbol}) =
    RecipeGraph(describe_graph(device, source))
# Layout stays separate from traversal: every input gets its own small block.
function graph_row!(lines, prefix, text, width)
    available = max(12, width - textwidth(prefix))
    words = split(text)
    row = ""
    for word in words
        if !isempty(row) && textwidth(row) + 1 + textwidth(word) > available
            push!(lines, prefix * row)
            row = String(word)
        else
            row = isempty(row) ? String(word) : row * " " * word
        end
    end
    return !isempty(row) && push!(lines, prefix * row)
end

# Render structured node metadata without parsing an encoded title.
function graph_node_rows!(lines, node, prefix, width; root = false)
    dimension = isnothing(node.dimension) ? "unknown" : string(node.dimension)
    if node.kind === :leaf
        root && graph_row!(lines, prefix, "target [leaf]", width)
        graph_row!(lines, prefix, "Leaf: " * string(node.component), width)
        graph_row!(lines, prefix, "Initial dimension: " * (isnothing(node.dimension) ? "default/unknown" : dimension), width)
    elseif node.kind === :combine
        name = isnothing(node.name) ? "" : " :" * string(node.name)
        heading = (root ? "target " : "") * "[combine]" * name
        graph_row!(lines, prefix, heading, width)
        graph_row!(lines, prefix, "Projected dimension: " * dimension, width)
    else
        graph_row!(lines, prefix, (root ? "target · " : "") * node.message, width)
    end
    for detail in node.details
        graph_row!(lines, prefix, detail, width)
    end
    return nothing
end

# Render ordered input edges with their selection summaries.
function graph_lines!(lines, node, prefix = "", width = 80)
    isempty(node.edges) && return lines
    push!(lines, rstrip(prefix))
    for (index, edge) in enumerate(node.edges)
        child = edge.child
        lastchild = index == length(node.edges)
        branchprefix = prefix * (lastchild ? "    " : "│   ")
        push!(lines, prefix * (lastchild ? "└── " : "├── ") * string(edge.name))
        graph_row!(lines, branchprefix, "Retained states: " * state_description(edge.states), width)
        graph_node_rows!(lines, child, branchprefix, width)
        graph_lines!(lines, child, branchprefix, width)
        !lastchild && push!(lines, rstrip(branchprefix))
    end
    return lines
end

# The graph and default recipe display share these rows and branch layout.
function Base.show(io::IO, ::MIME"text/plain", view::RecipeGraph)
    get(io, :compact, false) && return show(io, view)
    width = displaysize(io)[2]
    lines = String[]
    graph_row!(lines, "", "Recipe graph · planned build (root = target; branches = upstream inputs)", width)
    push!(lines, "")
    graph_node_rows!(lines, view.root, "", width; root = true)
    graph_lines!(lines, view.root, "", width)
    return print(io, join(lines, "\n"))
end

Base.show(io::IO, view::RecipeGraph) =
    print(io, "RecipeGraph(", view.root.kind, ", upstream inputs=", length(view.root.edges), ")")
