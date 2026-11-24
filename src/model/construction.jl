"""
    build(device, recipe; parameters=Dict(), dressed_threads=nothing, tracking_threads=missing)

Build a fresh construction-basis model. Operators, numerical Hamiltonians, and
exposed eigenstates share coordinates; diagonalization supplies spectra/labels
without rotating operators. Explicit selections create new bases from upstream
eigenvectors. Eigensystems remain eager except pure reductions, which inherit
them without a new solve. Explicit parameters override all recipe values.
`tracking_threads=missing` inherits each node’s effective dressed worker count;
an explicit positive Int overrides tracking independently.
`dressed_threads` overrides the sample worker count for every dressed node in
this build. It must be a positive `Int`, or `nothing` to honor each recipe's
`dressing_options.thread_count`. Workers are capped by available default-pool
threads and sample count. This execution setting does not alter recipe provenance.
"""
function build(
    device::QuantumDevice,
    recipe::Union{Recipe, Symbol};
    parameters = Dict{Num, Number}(),
    dressed_threads = nothing,
    tracking_threads = missing,
)
    validate_dressed_threads(dressed_threads)
    validate_tracking_threads(tracking_threads; inherit=true)
    node,_ = prepare_node(device,recipe,parameters)
    return assemble_node(node; dressed_threads, tracking_threads)
end
build(device::QuantumDevice, recipe::Union{Recipe, Symbol}, parameters; dressed_threads=nothing, tracking_threads=missing) =
    build(device, recipe; parameters, dressed_threads, tracking_threads)

"""
Return the original physical symbolic Hamiltonian without numerical construction.
Dependency and selection validation still applies; numerical `max_dimension` caps
are enforced by `build`, so oversized recipes can be inspected before refinement.
"""
function get_hamiltonian(device::QuantumDevice, recipe::Union{Recipe, Symbol})
    return resolve_node(device,recipe;enforce_limits=false).hamiltonian
end

# A leaf keeps generated operators in their original construction coordinates.
# Reference H is assembled from the same complete-product cache used downstream.
function assemble_node(node;dressed_threads=nothing,tracking_threads=missing)
    node.definition isa Leaf || return assemble_composite(node;dressed_threads,tracking_threads)
    recipe = node.definition
    generators,dimension = node.generators,node.dimension
    operators = Dict(ref=>generator(dimension) for (ref,generator) in generators)
    all(operator -> operator isa QuantumObject,values(operators)) ||
        throw(ArgumentError("generators must return QuantumObjects"))
    identity = qeye(dimension)
    all(operator -> size(operator) == (dimension,dimension) && operator.dimensions == identity.dimensions,
        values(operators)) || throw(DimensionMismatch(
            "generators must return square $dimension-state operators with matching QuantumToolbox dimensions"))
    products = leaf_products(operators,node.required_products,dimension)
    H = compile_hamiltonian(node.hamiltonian_terms,products,node.parameters,Set{Num}(),identity)
    energies,eigenvectors = diagonalize(H)
    eigensystem = leaf_eigensystem(only(node.owners),energies,eigenvectors)
    public = Dict{OperatorRef,QuantumObject}(ref=>products[((ref,false),)] for ref in node.public_operators)
    return NumericalModel(ResolvedRecipe(recipe,copy(node.parameters),Pair{Symbol,ResolvedRecipe}[]),
        H,eigensystem,public,node.hamiltonian,Set(node.owners),products)
end

function assemble_composite(node;dressed_threads=nothing,tracking_threads=missing)
    models = [assemble_node(child;dressed_threads,tracking_threads) for child in node.children]
    indices,projectors = selected_inputs(node,models)
    dimensions = length.(indices)
    identities = qeye.(dimensions)
    identity = ktensor(identities)
    products = projected_products(models,projectors,node.required_products)
    initial = 0*identity
    for j in eachindex(models)
        initial += ktensor([k == j ? spectrum_operator(models[j].spectrum[indices[j]]) : identities[k]
                            for k in eachindex(models)])
    end
    interaction = compile_hamiltonian(node.interaction_terms,products,node.parameters,Set{Num}(),identity;
        herm_check=0)
    H = hermiticity_cleanup(initial+interaction,DEFAULT_HERM_CHECK)
    eigensystem = composite_eigensystem(node.definition,models,indices,initial,interaction,H;
        dressed_threads,tracking_threads)
    public = Dict{OperatorRef,QuantumObject}(ref=>products[((ref,false),)] for ref in node.public_operators)
    tree = [alias=>model.recipe for (alias,model) in zip(node.aliases,models)]
    return NumericalModel(ResolvedRecipe(node.definition,copy(node.parameters),tree),
        H,eigensystem,public,node.hamiltonian,Set(node.owners),products)
end
