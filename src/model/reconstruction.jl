# Derived, read-only model preparation. No cache survives independently of its
# model, and rebuilding never mutates a cache shared with the original model.
struct _ProjectedModelTerm
    words
    factors
    operator
    coefficient
end

struct _ModelReconstruction
    dimensions
    local_hamiltonians
    local_eigensystems
    dressing
    local_expressions
    idle_terms
    interaction_terms
end

function _project_model_factor(component, dimension, word)
    parent = foldl((A, key) -> A * component.operators[key], word; init=qeye(component.dimension))
    return truncate(parent, projection(component.dimension, 1:dimension))
end

function _prepare_model_terms(expression, components, dimensions)
    sources = _operator_sources(components)
    return map(_symbolic_terms(expression)) do (word, coefficient)
        words = [Symbol[] for _ in components]
        for name in word
            haskey(sources, name) || throw(ArgumentError("Unknown operator: $name"))
            i, key = sources[name]
            push!(words[i], key)
        end
        factors = [_project_model_factor(c, d, w) for (c, d, w) in zip(components, dimensions, words)]
        _ProjectedModelTerm(words, factors, tensor(factors...), coefficient)
    end
end

function _refresh_model_terms(terms, components, dimensions, changed)
    return map(terms) do term
        affected = [i for i in changed if !isempty(term.words[i])]
        isempty(affected) && return term
        factors = copy(term.factors)
        for i in affected
            factors[i] = _project_model_factor(components[i], dimensions[i], term.words[i])
        end
        _ProjectedModelTerm(term.words, factors, tensor(factors...), term.coefficient)
    end
end

_compiled_model_terms(groups...) = [(t.operator, t.coefficient) for group in groups for t in group]

const _MODEL_INPUT_FIELDS = (:components, :interactions, :coupling_parameters,
    :max_dimension, :truncation_dimensions, :gates)

function _copy_model_fields(model, changes)
    fields = map(fieldnames(DeviceModel)) do key
        get(changes, key, getfield(model, key))
    end
    return DeviceModel(fields...)
end

function _rebuild_model(model::DeviceModel, changes)
    all(key -> key in _MODEL_INPUT_FIELDS, keys(changes)) ||
        throw(ArgumentError("Only model physical inputs and gates are writable"))
    changes = Dict(key => value for (key, value) in changes if value !== getfield(model, key))
    isempty(changes) && return model
    all(key -> key === :gates, keys(changes)) && return _copy_model_fields(model, changes)

    components = get(changes, :components, model.components)
    length(components) == length(model.components) &&
        all(a.name == b.name for (a, b) in zip(components, model.components)) ||
        throw(ArgumentError("Component updates must preserve names and order"))
    truncations = if haskey(changes, :truncation_dimensions)
        changes[:truncation_dimensions]
    elseif haskey(changes, :components)
        Dict(new => model.truncation_dimensions[old] for (old, new) in zip(model.components, components)
            if haskey(model.truncation_dimensions, old))
    else
        model.truncation_dimensions
    end
    interactions = get(changes, :interactions, model.interactions)
    coupling = get(changes, :coupling_parameters, model.coupling_parameters)
    max_dimension = get(changes, :max_dimension, model.max_dimension)
    gates = get(changes, :gates, model.gates)
    dimensions = _model_dimensions(components, truncations, max_dimension)
    cache = model.reconstruction
    structural = haskey(changes, :max_dimension) || haskey(changes, :truncation_dimensions) ||
        dimensions != cache.dimensions || any(a.dimension != b.dimension for (a,b) in zip(components, model.components))
    if structural
        result = make_model(components, interactions, coupling; max_dimension, truncation_dimensions=truncations)
        return _copy_model_fields(result, Dict(:gates => copy(gates)))
    end

    changed = findall(i -> components[i] !== model.components[i], eachindex(components))
    parameters = _model_parameters(components, coupling)
    local_expressions = copy(cache.local_expressions)
    H0s = copy(cache.local_hamiltonians)
    eigensystems = copy(cache.local_eigensystems)
    operators = Dict{Symbol,Any}(pairs(model.operators))
    sources = _operator_sources(components)
    for i in changed
        c, d = components[i], dimensions[i]
        local_expressions[i] = _promote_name_for_model(c.hamiltonian, Symbol(c.name))
        H0s[i] = truncate(numerical(c.hamiltonian, c.operators, c.parameters), projection(c.dimension, 1:d))
        eigensystems[i] = _local_eigensystem(H0s[i])
        # Remove obsolete operator names if a custom component changed its set.
        for name in keys(model.operators)
            startswith(String(name), c.name * "_") && !haskey(sources, name) && delete!(operators, name)
        end
        for (name, (owner, key)) in sources
            owner == i || continue
            operators[name] = wrap(truncate(c.operators[key], projection(c.dimension, 1:d)), i, dimensions)
        end
    end
    hamiltonian = sum(local_expressions) + interactions
    missing_parameters = setdiff(parameter_keys(hamiltonian), Set(keys(parameters)))
    isempty(missing_parameters) || throw(ArgumentError("Missing idle parameters: $missing_parameters"))
    any(v -> v isa Function, values(parameters)) && throw(ArgumentError("Idle model parameters must be static"))

    local_structure_changed = any(local_expressions[i] != cache.local_expressions[i] for i in changed)
    idle_terms = if isempty(changed)
        cache.idle_terms
    elseif local_structure_changed
        _prepare_model_terms(sum(local_expressions), components, dimensions)
    else
        _refresh_model_terms(cache.idle_terms, components, dimensions, changed)
    end
    interaction_terms = if interactions !== model.interactions
        _prepare_model_terms(interactions, components, dimensions)
    elseif isempty(changed)
        cache.interaction_terms
    else
        _refresh_model_terms(cache.interaction_terms, components, dimensions, changed)
    end
    terms = isempty(changed) && interactions === model.interactions ? model.compiled_terms :
        _compiled_model_terms(idle_terms, interaction_terms)
    dressing = isempty(changed) ? cache.dressing : _prepare_dressing(H0s, eigensystems)
    Hi = _evaluate_terms(_compiled_model_terms(interaction_terms), parameters)
    dressed = _track_dressing(dressing, Hi)
    reconstruction = _ModelReconstruction(dimensions, H0s, eigensystems, dressing,
        local_expressions, idle_terms, interaction_terms)
    return DeviceModel(components, interactions, coupling, max_dimension,
        isempty(changed) ? model.operators : (; operators...), parameters, hamiltonian,
        to_sparse(_evaluate_terms(terms, parameters)), dressed, copy(gates),
        Dict(truncations), terms, reconstruction)
end
