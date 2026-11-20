# Preserve each input's ordered factors, including repeated operators/adjoints.
local_product(word,owners) = Tuple(item for item in word if item[1].owner in owners)

# Collect complete products before numerical assembly; empty child factors need
# only retained identities and do not require upstream identity cache entries.
function collect_required_products!(node,downstream=Set{Tuple}())
    union!(node.required_products,downstream)
    union!(node.required_products,(word for (word,_) in node.hamiltonian_terms))
    union!(node.required_products,(((ref,false),) for ref in node.public_operators))
    for child in node.children
        localwords = Set{Tuple}()
        for word in node.required_products
            factor = local_product(word,child.owners)
            isempty(factor) || push!(localwords,factor)
        end
        collect_required_products!(child,localwords)
    end
    return node.required_products
end

# Keep integer requests zero-based; occupation tuples use the tracking owner's
# existing lookup rather than constructing a second label dictionary.
function resolve_selection(model,request,alias)
    isnothing(request) && return collect(eachindex(model.spectrum))
    first(request) isa Tuple || return request .+ 1
    tracking = model.eigensystem.tracking
    tracking.tracked || throw(ArgumentError("occupation selection for input $alias requires upstream tracking; enable dressing"))
    width = length(model.eigensystem.components)
    all(label -> length(label)==width,request) || throw(ArgumentError(
        "occupation selection for input $alias requires $width levels in component order $(model.eigensystem.components)"))
    lookup = getfield(tracking,:_lookup)
    missing = [label for label in request if !haskey(lookup,label)]
    isempty(missing) || throw(ArgumentError("occupation labels $missing are unavailable in input $alias; retain more upstream states"))
    return [tracking_row(tracking,label) for label in request]
end

# Selected physical eigenstates are projector columns. Keep input tensor dims
# and a flat retained output space, without dense rank-one outer products.
function projector_from_states(states)
    reference = first(states)
    all(ket -> ket.dimensions==reference.dimensions,states) ||
        throw(DimensionMismatch("projector ket dimensions differ"))
    columns = reduce(hcat,[ket.data for ket in states])
    return QuantumObject(columns;type=Operator(),dims=(Tuple(reference.dims[1]),(length(states),)))
end

# Ordinary tensor inputs retain request order. A pure reduction establishes its
# energy-ordered retained basis here, once, before H/projector/eigensystem consume it.
function selected_inputs(node,models)
    indices = [resolve_selection(model,request,alias)
               for (model,request,alias) in zip(models,node.selections,node.aliases)]
    if recipe_stage_operation(node.definition,[model.recipe.definition for model in models]) === :projection
        selected = only(indices)
        permutation = sortperm(only(models).spectrum[selected];alg=Base.Sort.MergeSort)
        indices[1] = selected[permutation]
    end
    projectors = [projector_from_states([model.states[i] for i in selected])
                  for (model,selected) in zip(models,indices)]
    return indices,projectors
end

# Multiply the complete ordered word in leaf coordinates before any projection.
function leaf_products(operators,required,dimension)
    products = Dict{Tuple,QuantumObject}()
    for word in required
        product = qeye(dimension)
        for (ref,adjoint_ref) in word
            operator = operators[ref]
            product *= adjoint_ref ? operator' : operator
        end
        products[word] = product
    end
    return products
end

# Partition full words once, project unique nonempty local words once, then reuse
# those factors for every tensor product. Empty factors are retained identities.
function projected_products(models,projectors,required)
    partitions = Dict{Tuple,Vector{Tuple}}()
    localwords = [Set{Tuple}() for _ in models]
    for word in required
        parts = Tuple[local_product(word,model.owners) for model in models]
        partitions[word] = parts
        for (j,part) in enumerate(parts)
            push!(localwords[j],part)
        end
    end
    factors = [Dict{Tuple,QuantumObject}() for _ in models]
    for (j,(model,projector)) in enumerate(zip(models,projectors))
        for word in localwords[j]
            factors[j][word] = isempty(word) ? qeye(size(projector,2)) :
                projector'*model.operator_products[word]*projector
        end
    end
    products = Dict{Tuple,QuantumObject}()
    for (word,parts) in partitions
        products[word] = ktensor([factors[j][part] for (j,part) in enumerate(parts)])
    end
    return products
end
