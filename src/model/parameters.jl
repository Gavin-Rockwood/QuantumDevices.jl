# Numerical builds attach the same symbolic assignments used by inventory inspection.
function prepare_parameters!(node, parameters; components::AbstractDict, extra_terms=Interaction[], expressions=Sym[])
    resolution = resolve_recipe_inputs(node.definition,components;overrides=parameters,extra_terms,expressions)
    require_resolved_inputs(resolution)
    assign_node_inputs!(node,Dict(entry.path=>entry.values for entry in resolution.nodes))
    return node,Dict{Num,Number}(resolution.values)
end
function assign_node_inputs!(node, assignments; path=())
    node.parameters = Dict{Num,Number}(assignments[path])
    for (alias,child) in zip(node.aliases,node.children)
        assign_node_inputs!(child,assignments;path=(path...,alias))
    end
    return nothing
end
