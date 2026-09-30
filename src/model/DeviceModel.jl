"""
    DeviceModel

Coupled device built by [`make_model`](@ref). Exposes `components`,
`coupling_parameters`, symbolic `hamiltonian`, embedded `operators`, merged
`parameters`, sparse idle Hamiltonian `H`, `truncation_dimensions`, and `gates`.
`states`, `others`, and `confidence` forward to its dressed-state tracking result.

Use [`make_model`](@ref) rather than the positional field constructor. Derived
matrices and eigensystems are reconstructed by supported [`setpath`](@ref) updates.
`model.gates` is a mutable dictionary accepting named [`DeviceGate`](@ref) objects.
"""
struct DeviceModel
    components
    interactions
    coupling_parameters
    max_dimension
    operators
    parameters
    hamiltonian
    H
    eigensystem
    gates
    truncation_dimensions
    compiled_terms
end

function _promote_name_for_model(x, name::Symbol)
    x isa Sym || return x

    e = x.expr

    if e[1] === :op
        return op(Symbol(name, "_", e[2]))

    elseif e[1] === :param
        return param(Symbol(name, "_", e[2]))

    elseif e[1] === :val
        return x

    elseif e[1] === :call
        return call(
            e[2],
            map(a -> _promote_name_for_model(a, name), e[3:end])...
        )
    end

    error("Unknown expression: $e")
end

function _model_dimensions(components, truncation_dimensions, max_dimension)
    isempty(components) && throw(ArgumentError("A model needs at least one component"))
    allunique(c.name for c in components) || throw(ArgumentError("Component names must be unique"))
    all(c -> c in components, keys(truncation_dimensions)) || throw(ArgumentError("Truncation refers to an unknown component"))
    dims = Int[]
    for component in components
        component.dimension isa Integer && component.dimension > 0 || throw(ArgumentError("Invalid parent dimension"))
        d = get(truncation_dimensions, component, component.dimension)
        d isa Integer && 1 <= d <= component.dimension || throw(ArgumentError("Invalid retained dimension for $(component.name)"))
        push!(dims, d)
    end
    max_dimension > 0 || throw(ArgumentError("max_dimension must be positive"))
    prod(big.(dims)) <= max_dimension || throw(ArgumentError("Total dimension exceeds max_dimension; provide truncation_dimensions or raise the limit"))
    return dims
end

function _operator_sources(components)
    sources = Dict{Symbol,Tuple{Int,Symbol}}()
    for (i, component) in enumerate(components), key in keys(component.operators)
        name = Symbol(component.name, "_", key)
        haskey(sources, name) && throw(ArgumentError("Duplicate operator name: $name"))
        size(component.operators[key]) == (component.dimension, component.dimension) || throw(DimensionMismatch("Invalid operator size for $name"))
        sources[name] = (i, key)
    end
    return sources
end

function _projected_terms(expression, components, dims)
    sources = _operator_sources(components)
    projectors = [projection(c.dimension, 1:d) for (c, d) in zip(components, dims)]
    return map(_symbolic_terms(expression)) do (word, coefficient)
        local_words = [Symbol[] for _ in components]
        for name in word
            haskey(sources, name) || throw(ArgumentError("Unknown operator: $name"))
            i, key = sources[name]
            push!(local_words[i], key)
        end
        factors = map(eachindex(components)) do i
            c = components[i]
            # Multiply in the parent basis, including repeated operators, before truncation.
            parent = foldl((A, key) -> A * c.operators[key], local_words[i]; init = qeye(c.dimension))
            truncate(parent, projectors[i])
        end
        (tensor(factors...), coefficient)
    end
end

"""
    make_model(components, interaction::Sym, coupling_parameters::NamedTuple;
               max_dimension=10^4, truncation_dimensions=Dict())

Build a coupled model in the supplied tensor-product component order. Names must
be unique. Component parameters/operators become `componentname_localname`.
Coupling names must not collide with them; all idle parameters must be static.

`truncation_dimensions` maps component objects to retained dimensions between one
and their parent dimensions. Their product must not exceed `max_dimension`.
Complete ordered operator products are multiplied in the parent basis before
projection, preserving virtual excursions omitted by multiplying truncated factors.

Returns a [`DeviceModel`](@ref) with a sparse idle Hamiltonian and dressed states
tracked from the uncoupled system. Construction rejects unknown operators,
unresolved parameters, invalid truncations, and naming collisions.
"""
function make_model(components, interaction::Sym, coupling_parameters::NamedTuple;
                    max_dimension = 10^4, truncation_dimensions = Dict())
    components = collect(components)
    dims = _model_dimensions(components, truncation_dimensions, max_dimension)
    sources = _operator_sources(components)
    parameters = Dict{Symbol,Any}()
    for component in components, key in keys(component.parameters)
        name = Symbol(component.name, "_", key)
        haskey(parameters, name) && throw(ArgumentError("Duplicate parameter name: $name"))
        parameters[name] = component.parameters[key]
    end
    for (key, value) in pairs(coupling_parameters)
        haskey(parameters, key) && throw(ArgumentError("Coupling parameter collides with component parameter: $key"))
        parameters[key] = value
    end
    parameters = (; parameters...)
    operators = Dict{Symbol,Any}()
    for (name, (i, key)) in sources
        c = components[i]
        operators[name] = wrap(truncate(c.operators[key], projection(c.dimension, 1:dims[i])), i, dims)
    end
    operators = (; operators...)

    local_hamiltonians = [_promote_name_for_model(c.hamiltonian, Symbol(c.name)) for c in components]

    hamiltonian = sum(local_hamiltonians) + interaction

    missing_parameters = setdiff(parameter_keys(hamiltonian), Set(keys(parameters)))
    isempty(missing_parameters) || throw(ArgumentError("Missing idle parameters: $missing_parameters"))
    any(v -> v isa Function, values(parameters)) && throw(ArgumentError("Idle model parameters must be static"))

    terms = _projected_terms(hamiltonian, components, dims)

    H0s = [truncate(numerical(c.hamiltonian, c.operators, c.parameters), projection(c.dimension, 1:d)) for (c, d) in zip(components, dims)]
    Hi = _evaluate_terms(_projected_terms(interaction, components, dims), parameters)
    dressed = get_dressed_states(H0s, Hi)

    return DeviceModel(components, interaction, coupling_parameters, max_dimension,
        operators, parameters, hamiltonian,
        to_sparse(_evaluate_terms(terms, parameters)), dressed, Dict{Any,Any}(),
        Dict(truncation_dimensions), terms)
end
