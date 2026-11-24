# Capture resolves names and validates structure; this traversal only prepares
# numerical configuration, dimensions, stable requests, and physical terms.
function resolve_node(device,source;enforce_limits::Bool=true)
    components = Dict{Symbol,Component}()
    captured = snapshot_recipe(device,source;components)
    return prepare_captured_node(captured,components;enforce_limits)
end

physical_terms(expression) = Tuple{Tuple,Sym}[(Tuple(word),coefficient) for (word,coefficient) in terms(expression)]

function prepare_captured_node(definition,components;enforce_limits=true,path="target")
    if definition isa Leaf
        component = components[definition.component]
        generators = merge(component_operator_generators(component),definition.generators)
        dimension = isnothing(definition.dimension) ? component.dimension : definition.dimension
        isnothing(dimension) && throw(ArgumentError("$path: specify Leaf.dimension or Component.dimension"))
        dimension = validate_operator_dimension(dimension)
        return BuildNode(;definition,owners=[component.name],generators,
            relevant=Set{Num}(values(component.variables)),dimension,
            hamiltonian=component.hamiltonian,hamiltonian_terms=physical_terms(component.hamiltonian),
            public_operators=collect(OperatorRef,keys(generators)))
    end
    inputs = recipe_inputs(definition)
    aliases = [input.name for input in inputs]
    sources = [input.recipe for input in inputs]
    selections = Any[input.states for input in inputs]
    children = [prepare_captured_node(source,components;enforce_limits,path="$path -> $alias")
                for (alias,source) in zip(aliases,sources)]
    limit = definition.max_dimension
    limit >= 0 || throw(ArgumentError("$path: max_dimension must be nonnegative"))
    dimension = prod((big(isnothing(request) ? child.dimension : length(request))
                      for (child,request) in zip(children,selections));init=big(1))
    enforce_limits && limit>0 && dimension>limit &&
        throw(ArgumentError("$path: requested composite dimension $dimension exceeds max_dimension $limit"))
    dimension <= typemax(Int) || throw(ArgumentError(
        "$path: requested composite dimension $dimension exceeds representable dimension $(typemax(Int))"))
    for (alias,child,request) in zip(aliases,children,selections)
        isnothing(request) || first(request) isa Tuple || all(i->i<child.dimension,request) ||
            throw(ArgumentError("$path -> $alias: retained state out of bounds for input dimension $(child.dimension)"))
    end
    interactions = collect(Interaction,definition.interactions)
    physical_interactions = physical_expression.([interaction.expression for interaction in interactions])
    interaction_terms = reduce(vcat,physical_terms.(physical_interactions);init=Tuple{Tuple,Sym}[])
    hamiltonian_terms = reduce(vcat,(child.hamiltonian_terms for child in children);init=Tuple{Tuple,Sym}[])
    append!(hamiltonian_terms,interaction_terms)
    hamiltonian = foldl(+,physical_interactions;init=foldl(+,(child.hamiltonian for child in children)))
    relevant = union((child.relevant for child in children)...,
        Set{Num}(variable for interaction in interactions for variable in values(interaction.variables)))
    return BuildNode(;definition,aliases,children,selections,interactions,relevant,
        owners=reduce(vcat,(child.owners for child in children)),dimension=Int(dimension),
        hamiltonian,hamiltonian_terms,interaction_terms,
        public_operators=reduce(vcat,(child.public_operators for child in children);init=OperatorRef[]))
end

# Prepare a captured reference, its operating point, and required model products.
function prepare_node(device,recipe,parameters)
    validate_component_input_defaults(device.components)
    components = Dict{Symbol,Component}()
    captured = snapshot_recipe(device,recipe;components)
    node = prepare_captured_node(captured,components)
    node,resolved = prepare_parameters!(node,parameters;components)
    collect_required_products!(node)
    return node,resolved
end
