"""
    DeviceModel

Coupled device built by [`make_model`](@ref). Exposes `components`,
`coupling_parameters`, symbolic `hamiltonian`, embedded `operators`, merged
`parameters`, sparse idle Hamiltonian `H`, `truncation_dimensions`, and `gates`.
Hamiltonians and dimensional energy parameters use frequency units (`E/h`,
cycles per unit time); `numerical` does not multiply them by `2π`.
`states`, `others`, and `confidence` forward to its dressed-state tracking result.
`spectrum` is the label-to-energy dictionary, sharing the same keys as `states`;
energies are in frequency units, `E/h` (cycles per unit time).

Use [`make_model`](@ref) rather than the positional field constructor. Derived
matrices and eigensystems are reconstructed by supported [`setpath`](@ref) updates.
Gate updates preserve idle physics. Physical updates use the normal builder and
do not inherit gates. Reference tracking can reuse dressed-state labels.
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
               max_dimension=10^4, truncation_dimensions=Dict(),
               reference_model=nothing, tracking_steps=nothing, threaded_components=true)

Build a coupled model in the supplied tensor-product component order. Names must
be unique. Component parameters/operators become `componentname_localname`.
Coupling names must not collide with them; all idle parameters must be static.

`truncation_dimensions` maps component objects to retained dimensions between one
and their parent dimensions. Their product must not exceed `max_dimension`.
Complete ordered operator products are multiplied in the parent basis before
projection, preserving virtual excursions omitted by multiplying truncated factors.

Returns a [`DeviceModel`](@ref) with a sparse idle Hamiltonian and dressed states
tracked from the uncoupled system, or from `reference_model`. `tracking_steps`
counts intermediate models: zero diagonalizes the final Hamiltonian once and
matches its eigenbasis directly. Defaults are zero with a reference and 18
without one (preserving the original 20-point bare-state trajectory).
Reference tracking transports component bases to compare physical overlaps.
Positive steps interpolate numeric inputs and require unchanged model structure.
Models start with an empty gate dictionary. `threaded_components=true` prepares
local Hamiltonians concurrently on Julia threads. Construction rejects unknown operators,
unresolved parameters, invalid truncations, and naming collisions.
"""
function make_model(components, interaction::Sym, coupling_parameters::NamedTuple;
                    max_dimension = 10^4, truncation_dimensions = Dict(),
                    reference_model=nothing,tracking_steps=nothing,threaded_components::Bool=true)
    components = collect(components)
    dims = _model_dimensions(components, truncation_dimensions, max_dimension)
    sources = _operator_sources(components)
    parameters = _model_parameters(components, coupling_parameters)
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

    terms = _projected_terms(hamiltonian,components,dims)
    H = to_sparse(_evaluate_terms(terms,parameters))
    steps = _tracking_step_count(tracking_steps,reference_model)
    dressed = if reference_model === nothing
        H0s = Vector{Any}(undef,length(components))
        function prepare_local(i)
            c,d = components[i],dims[i]
            H0s[i] = truncate(numerical(c.hamiltonian,c.operators,c.parameters),projection(c.dimension,1:d))
        end
        if threaded_components && length(components)>1 && Threads.nthreads()>1
            Threads.@threads for i in eachindex(components)
                prepare_local(i)
            end
        else
            foreach(prepare_local,eachindex(components))
        end
        preparation = _prepare_dressing(H0s)
        _track_dressing(preparation,H-preparation.H0;nsteps=steps+2)
    else
        _track_model_reference(components,interaction,coupling_parameters,H,dims,reference_model,steps;
            max_dimension,threaded_components)
    end
    return DeviceModel(components,interaction,coupling_parameters,max_dimension,
        operators,parameters,hamiltonian,H,dressed,Dict{Any,Any}(),Dict(truncation_dimensions),terms)
end

function _model_parameters(components, coupling_parameters)
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
    return parameters
end
