"""
    InteractionSpec(id, expr; metadata = Dict{Symbol,Any}())

Symbolic interaction term between components in a `ModelSpec`.
"""
struct InteractionSpec
    id::Symbol
    expr::OperatorExpr
    metadata::Dict{Symbol,Any}
end

function InteractionSpec(
        id::Symbol,
        expr::OperatorExpr;
        metadata = Dict{Symbol,Any}())
    return InteractionSpec(id, expr, _metadata_dict(metadata))
end

"""
    DressingSpec(; steps = 20, schedule = :power, exponent = 3.0,
                  minimum_overlap = 0.5, on_low_overlap = :warn)

Options for adiabatically tracking product states into interacting model states.
"""
struct DressingSpec
    steps::Int
    schedule::Symbol
    exponent::Float64
    minimum_overlap::Float64
    on_low_overlap::Symbol
end

function DressingSpec(;
        steps::Integer = 20,
        schedule::Symbol = :power,
        exponent::Real = 3.0,
        minimum_overlap::Real = 0.5,
        on_low_overlap::Symbol = :warn)
    steps >= 2 || error("DressingSpec steps must be at least two.")
    schedule in (:linear, :power) ||
        error("Unknown dressing schedule $schedule.")
    exponent > 0 || error("DressingSpec exponent must be positive.")
    0 <= minimum_overlap <= 1 ||
        error("DressingSpec minimum_overlap must be between zero and one.")
    on_low_overlap in (:warn, :error, :ignore) ||
        error("Unknown DressingSpec low-overlap policy $on_low_overlap.")
    return DressingSpec(
        Int(steps),
        schedule,
        Float64(exponent),
        Float64(minimum_overlap),
        on_low_overlap,
    )
end

"""
    ModelSpec(id; components, children, interactions, dims, ...)

Declarative model tree with local components, child models, interactions,
truncation dimensions, and optional selected states.
"""
struct ModelSpec
    id::Symbol
    components::Dict{Symbol,Component}
    children::Dict{Symbol,ModelSpec}
    interactions::Dict{Symbol,InteractionSpec}
    dims::Dict{Symbol,Int}
    initialization_dims::Dict{Any,Int}
    states_to_keep::Vector{Tuple{Vararg{Int}}}
    dressing::DressingSpec
    parameters::Dict{Symbol,DeviceParameter}
    metadata::Dict{Symbol,Any}
end

function ModelSpec(
        id::Symbol;
        components = Dict{Symbol,Component}(),
        children = Dict{Symbol,ModelSpec}(),
        interactions = Dict{Symbol,InteractionSpec}(),
        dims,
        initialization_dims = Dict{Any,Int}(),
        states_to_keep = Tuple[],
        dressing::DressingSpec = DressingSpec(),
        parameters = Dict{Symbol,DeviceParameter}(),
        metadata = Dict{Symbol,Any}())
    component_dict = Dict{Symbol,Component}(components)
    child_dict = Dict{Symbol,ModelSpec}(children)
    interaction_dict = Dict{Symbol,InteractionSpec}(interactions)
    dimension_dict = Dict{Symbol,Int}(dims)
    initialization_dict = Dict{Any,Int}(initialization_dims)
    selected_states = Tuple{Vararg{Int}}[
        Tuple(Int.(collect(label)))
        for label in states_to_keep
    ]

    overlap = intersect(keys(component_dict), keys(child_dict))
    isempty(overlap) ||
        error("ModelSpec $id has duplicate component and child IDs: $(join(sort!(collect(overlap)), ", ")).")

    available = union(Set(keys(component_dict)), Set(keys(child_dict)))
    Set(keys(dimension_dict)) == available ||
        error("ModelSpec $id dims must contain exactly its component and child IDs.")
    all(dimension -> dimension > 0, values(dimension_dict)) ||
        error("All dimensions for ModelSpec $id must be positive.")
    length(unique(selected_states)) == length(selected_states) ||
        error("ModelSpec $id states_to_keep contains duplicate tuples.")
    subsystem_order = sort!(collect(available))
    for label in selected_states
        length(label) == length(subsystem_order) ||
            error("State tuple $label for ModelSpec $id must have $(length(subsystem_order)) entries.")
        for (index, subsystem_id) in enumerate(subsystem_order)
            0 <= label[index] < dimension_dict[subsystem_id] ||
                error("State tuple $label index $(label[index]) is outside subsystem $subsystem_id dimension $(dimension_dict[subsystem_id]).")
        end
    end

    for (component_id, component) in component_dict
        component.id == component_id ||
            error("Component dictionary key $component_id does not match component ID $(component.id).")
    end
    for (child_id, child) in child_dict
        child.id == child_id ||
            error("Child dictionary key $child_id does not match child ModelSpec ID $(child.id).")
    end
    for (path, dimension) in initialization_dict
        dimension > 0 || error("Initialization dimension for $path must be positive.")
        _validate_initialization_path(id, component_dict, child_dict, path)
    end
    for interaction in values(interaction_dict)
        participants = _participating_components(interaction.expr)
        missing = setdiff(participants, available)
        isempty(missing) ||
            error("Interaction $(interaction.id) references components not present in ModelSpec $id: $(join(sort!(collect(missing)), ", ")).")
    end

    return ModelSpec(
        id,
        component_dict,
        child_dict,
        interaction_dict,
        dimension_dict,
        initialization_dict,
        selected_states,
        dressing,
        Dict{Symbol,DeviceParameter}(parameters),
        _metadata_dict(metadata),
    )
end

"""
    _validate_initialization_path(model_id, components, children, path)

Validate that an initialization-dimension override points to an existing local
component or nested child component.
"""
function _validate_initialization_path(model_id, components, children, path)
    if path isa Symbol
        (haskey(components, path) || haskey(children, path)) ||
            error("Initialization path $path does not refer to a component or generated child of ModelSpec $model_id.")
        return
    end
    path isa Tuple{Vararg{Symbol}} && !isempty(path) ||
        error("Initialization paths must be symbols or nonempty symbol tuples.")
    child_id = path[1]
    haskey(children, child_id) ||
        error("Initialization path $(join(path, "/")) does not begin with a child of ModelSpec $model_id.")
    length(path) > 1 ||
        error("Nested initialization path $(join(path, "/")) must identify a component inside child $child_id.")
    child = children[child_id]
    tail = length(path) == 2 ? path[2] : path[2:end]
    _validate_initialization_path(child.id, child.components, child.children, tail)
end

"""
    QuantumDeviceModel

Numerical model produced from a `ModelSpec`, including Hamiltonian, embedded
operators, dressed/product states, energies, and provenance.
"""
struct QuantumDeviceModel
    hamiltonian
    operators::Dict{Tuple{Vararg{Symbol}},Any}
    states
    energies
    resolved_parameters::Dict{Symbol,Any}
    spec::ModelSpec
    basis::Symbol
    state_order::Tuple{Vararg{Symbol}}
    state_labels::Vector{Tuple{Vararg{Int}}}
    state_overlaps::Matrix{Float64}
end

"""
    FrozenModelSpec

Component spec produced by freezing a `QuantumDeviceModel` into a reusable
component.
"""
struct FrozenModelSpec <: AbstractComponentSpec
    dimension::DeviceParameter
    operators::Dict{Any,Any}
    hamiltonian::OperatorExpr
end

"""
    numerical_operator(spec::FrozenModelSpec, operator; dimension = nothing)

Return an already-projected operator stored on a frozen model component.
"""
function numerical_operator(
        spec::FrozenModelSpec,
        operator,
        args...;
        dimension = nothing,
        kwargs...)
    haskey(spec.operators, operator) ||
        error("Operator $operator not in frozen model operators.")
    expected = spec.dimension.default
    dimension === nothing || dimension == expected ||
        error("Frozen model dimension is $expected, received $dimension.")
    return spec.operators[operator]
end

"""
    model(spec::ModelSpec; params = Dict(), kwargs...)

Build a numerical `QuantumDeviceModel` from a symbolic model tree.
"""
model(spec::ModelSpec; params = Dict(), kwargs...) =
    _build_model(spec, Dict{Any,Int}(); params = params, kwargs...)

"""
    _build_model(spec, inherited_initialization; params, kwargs...)

Materialize children, embed local operators, assemble the Hamiltonian, and
track product states for one model-tree node.
"""
function _build_model(
        spec::ModelSpec,
        inherited_initialization::Dict{Any,Int};
        params,
        kwargs...)
    initialization = merge(copy(spec.initialization_dims), inherited_initialization)
    resolved_spec = _resolve_model_spec(spec, initialization)
    resolved_parameters = _resolve_model_parameters(resolved_spec, params, kwargs)
    components, resolved_children = _materialize_components(
        resolved_spec,
        initialization;
        params = params,
        kwargs...,
    )
    resolved_spec = ModelSpec(
        resolved_spec.id;
        components = resolved_spec.components,
        children = resolved_children,
        interactions = resolved_spec.interactions,
        dims = resolved_spec.dims,
        initialization_dims = resolved_spec.initialization_dims,
        states_to_keep = resolved_spec.states_to_keep,
        dressing = resolved_spec.dressing,
        parameters = resolved_spec.parameters,
        metadata = resolved_spec.metadata,
    )
    isempty(components) && error("ModelSpec $(resolved_spec.id) has no components.")

    order = sort!(collect(keys(components)))
    dimensions = resolved_spec.dims
    local_hamiltonians = Dict(
        id => numerical(
            components[id];
            dimension = dimensions[id],
            params = params,
            kwargs...,
        )
        for id in order
    )
    operators = _embedded_operators(
        components,
        order,
        dimensions;
        params = params,
        kwargs...,
    )
    uncoupled_hamiltonian = reduce(
        +,
        (
            _embed_operator(local_hamiltonians[id], id, order, dimensions)
            for id in order
        ),
    )
    interaction_hamiltonian = isempty(resolved_spec.interactions) ?
        0 * uncoupled_hamiltonian :
        reduce(
            +,
            (
                _numerical_model_expr(interaction.expr, operators, resolved_parameters)
                for interaction in values(resolved_spec.interactions)
            ),
        )
    hamiltonian = uncoupled_hamiltonian + interaction_hamiltonian

    labels = isempty(resolved_spec.states_to_keep) ?
        _all_state_labels(order, dimensions) :
        resolved_spec.states_to_keep
    references = _reference_product_states(
        components,
        local_hamiltonians,
        order,
        labels,
    )
    dressing = resolved_spec.dressing
    tracking = get_dressed_states(
        uncoupled_hamiltonian,
        interaction_hamiltonian,
        references;
        steps = dressing.steps,
        schedule = dressing.schedule,
        exponent = dressing.exponent,
        minimum_overlap = dressing.minimum_overlap,
        on_low_overlap = dressing.on_low_overlap,
    )
    dressed_states = Dict(
        label => QuantumObject(tracking.states[:, index])
        for (index, label) in enumerate(labels)
    )
    dressed_energies = Dict(
        label => real(tracking.energies[index])
        for (index, label) in enumerate(labels)
    )

    if isempty(resolved_spec.states_to_keep)
        return QuantumDeviceModel(
            hamiltonian,
            operators,
            dressed_states,
            dressed_energies,
            resolved_parameters,
            resolved_spec,
            :product,
            Tuple(order),
            copy(labels),
            tracking.overlaps,
        )
    end

    dressed_operators = Dict(
        path => _project_operator(tracking.states, operator)
        for (path, operator) in operators
    )
    dressed_hamiltonian = QuantumObject(
        Diagonal(ComplexF64.(tracking.energies)),
    )
    return QuantumDeviceModel(
        dressed_hamiltonian,
        dressed_operators,
        dressed_states,
        dressed_energies,
        resolved_parameters,
        resolved_spec,
        :dressed,
        Tuple(order),
        copy(labels),
        tracking.overlaps,
    )
end

"""
    _resolve_model_spec(spec, initialization)

Turn oversized or infinite local components into child model specs so each node
builds from a valid initialization dimension.
"""
function _resolve_model_spec(spec::ModelSpec, initialization)
    components = copy(spec.components)
    children = copy(spec.children)

    for id in sort!(collect(keys(spec.components)))
        component = spec.components[id]
        requested = spec.dims[id]
        haskey(component.parameters, :dimension) ||
            error("Component $(component.id) does not define a maximum dimension.")
        maximum = component.parameters[:dimension].default
        initial = if isfinite(maximum)
            Int(maximum)
        else
            haskey(initialization, id) ||
                error("Component $id has infinite maximum dimension; define initialization_dims[$(repr(id))].")
            initialization[id]
        end

        requested <= initial ||
            error("ModelSpec $(spec.id) requests dimension $requested for $id from initialization dimension $initial.")
        if component.spec isa Union{TransmonSpec,FluxTunableTransmonSpec}
            isodd(initial) ||
                error("Transmon-like initialization dimension for $id must be odd.")
        end

        if requested < initial
            pop!(components, id)
            children[id] = ModelSpec(
                id;
                components = Dict(id => component),
                dims = Dict(id => initial),
                initialization_dims = isfinite(maximum) ?
                    Dict{Any,Int}() :
                    Dict{Any,Int}(id => initial),
            )
        end
    end

    return ModelSpec(
        spec.id;
        components = components,
        children = children,
        interactions = spec.interactions,
        dims = spec.dims,
        initialization_dims = spec.initialization_dims,
        states_to_keep = spec.states_to_keep,
        dressing = spec.dressing,
        parameters = spec.parameters,
        metadata = spec.metadata,
    )
end

"""
    component(model::QuantumDeviceModel)

Project a numerical model into a frozen component with the model's full
dimension.
"""
function component(model::QuantumDeviceModel)
    return _component_from_model(model, size(model.hamiltonian, 1))
end

"""
    _component_from_model(model, dimension)

Project a model's Hamiltonian and operators into a frozen component basis of
the requested dimension.
"""
function _component_from_model(model::QuantumDeviceModel, dimension::Int)
    full_dimension = size(model.hamiltonian, 1)
    1 <= dimension <= full_dimension ||
        error("Projected model dimension must be between 1 and $full_dimension.")
    projection = model.basis == :product ?
        _model_state_matrix(model, dimension) :
        Matrix{ComplexF64}(I, full_dimension, dimension)
    projected = Dict{Any,Any}()
    single_component = length(model.spec.components) + length(model.spec.children) == 1
    for (path, operator) in model.operators
        public_path = single_component ? path[2:end] : path
        key = length(public_path) == 1 ? public_path[1] : public_path
        projected[key] = _project_operator(projection, operator)
    end
    projected[:hamiltonian] = _project_operator(projection, model.hamiltonian)

    dimension_parameter = DeviceParameter(
        ParamPath(:dimension);
        domain = integerrange(dimension),
        fixed = true,
        required = true,
        default = dimension,
        description = "Maximum dressed-model dimension",
    )
    frozen = FrozenModelSpec(dimension_parameter, projected, op(:hamiltonian))
    return Component(
        frozen,
        model.spec.id,
        _new_uuid(),
        String(model.spec.id),
        merge(
            copy(model.spec.metadata),
            Dict{Symbol,Any}(
                :source_model => model.spec.id,
                :resolved_parameters => copy(model.resolved_parameters),
            ),
        ),
    )
end

"""
    _model_state_matrix(model, dimension)

Construct the projection matrix from the lowest-energy product-basis states in
a product-basis model.
"""
function _model_state_matrix(model::QuantumDeviceModel, dimension)
    labels = sort(
        model.state_labels;
        by = label -> (
            model.energies[label],
            findfirst(==(label), model.state_labels),
        ),
    )
    selected = labels[1:dimension]
    return hcat((vec(model.states[label].data) for label in selected)...)
end

"""
    _materialize_components(spec, initialization; params, kwargs...)

Build child models, freeze them back into components, and return the component
catalog used by this model node.
"""
function _materialize_components(spec, initialization; params, kwargs...)
    components = copy(spec.components)
    resolved_children = Dict{Symbol,ModelSpec}()
    for (id, child) in spec.children
        child_initialization = Dict{Any,Int}()
        for (path, dimension) in initialization
            path isa Tuple || continue
            !isempty(path) && path[1] == id || continue
            tail = length(path) == 2 ? path[2] : path[2:end]
            child_initialization[tail] = dimension
        end
        child_model = _build_model(
            child,
            child_initialization;
            params = params,
            kwargs...,
        )
        child_dimension = size(child_model.hamiltonian, 1)
        requested = spec.dims[id]
        requested <= child_dimension ||
            error("ModelSpec $(spec.id) requests dimension $requested from child $id with dimension $child_dimension.")
        components[id] = _component_from_model(child_model, requested)
        resolved_children[id] = child_model.spec
    end
    return components, resolved_children
end

"""
    _resolve_model_parameters(spec, params, kwargs)

Resolve model-level parameters from keyword overrides, `params`, or defaults.
"""
function _resolve_model_parameters(spec::ModelSpec, params, kwargs)
    values = Dict{Symbol,Any}()
    explicit = Dict{Symbol,Any}(pairs(kwargs))
    for (key, parameter) in spec.parameters
        value = if haskey(explicit, key)
            explicit[key]
        elseif haskey(params, key)
            params[key]
        elseif haskey(params, parameter.path)
            params[parameter.path]
        else
            parameter.default
        end
        value in parameter.domain ||
            error("Value $value is outside the domain for model parameter $key.")
        values[key] = value
    end
    return values
end

"""
    _embedded_operators(components, order, dimensions; params, kwargs...)

Evaluate each component operator and embed it into the full tensor-product
model space.
"""
function _embedded_operators(components, order, dimensions; params, kwargs...)
    operators = Dict{Tuple{Vararg{Symbol}},Any}()
    for id in order
        component = components[id]
        for operator in component.operators
            local_path = operator isa Symbol ? (operator,) : operator
            value = numerical(
                component,
                operator;
                dimension = dimensions[id],
                params = params,
                kwargs...,
            )
            operators[(id, local_path...)] =
                _embed_operator(value, id, order, dimensions)
        end
    end
    return operators
end

"""
    _reference_product_states(components, local_hamiltonians, order, labels)

Build bare product states used as references for dressed-state tracking.
"""
function _reference_product_states(components, local_hamiltonians, order, labels)
    local_states = Dict(
        id => components[id].spec isa FrozenModelSpec ?
            Matrix{ComplexF64}(I, size(local_hamiltonians[id], 1), size(local_hamiltonians[id], 1)) :
            eigenstates(local_hamiltonians[id]).vectors
        for id in order
    )
    product_dimension = prod(size(local_states[id], 1) for id in order)
    references = Matrix{ComplexF64}(undef, product_dimension, length(labels))
    for (column, label) in enumerate(labels)
        state = ComplexF64[1]
        for (index, id) in enumerate(order)
            state = kron(state, local_states[id][:, label[index] + 1])
        end
        references[:, column] = state
    end
    return references
end

"""
    _all_state_labels(order, dimensions)

Enumerate product-basis state labels in the model's subsystem order.
"""
function _all_state_labels(order, dimensions)
    ranges = (0:(dimensions[id] - 1) for id in order)
    labels = Tuple{Vararg{Int}}[
        Tuple(label)
        for label in Iterators.product(ranges...)
    ]
    return vec(labels)
end

"""
    _numerical_model_expr(expr, operators, parameters)

Evaluate a model-level operator expression against embedded operators and
resolved model parameters.
"""
function _numerical_model_expr(expr::OperatorExpr, operators, parameters)
    if expr.head == :op
        path = expr.args[1].path
        length(path) >= 4 && path[1] == :components && path[3] == :operators ||
            error("Model operator paths must be rooted at a component.")
        key = (path[2], path[4:end]...)
        haskey(operators, key) ||
            error("Operator path $(join(key, "/")) is not available in the model.")
        return operators[key]
    elseif expr.head == :+
        return reduce(+, (_numerical_model_expr(arg, operators, parameters) for arg in expr.args))
    elseif expr.head == :*
        return reduce(*, (_numerical_model_factor(arg, operators, parameters) for arg in expr.args))
    elseif expr.head == :^
        return _numerical_model_expr(expr.args[1], operators, parameters) ^ expr.args[2]
    elseif expr.head == :adjoint
        return adjoint(_numerical_model_expr(expr.args[1], operators, parameters))
    end
    error("Unknown operator expression head $(expr.head).")
end

_numerical_model_factor(expr::OperatorExpr, operators, parameters) =
    _numerical_model_expr(expr, operators, parameters)
_numerical_model_factor(expr::ScalarExpr, operators, parameters) =
    _numerical_model_scalar(expr, parameters)
_numerical_model_factor(value::Number, operators, parameters) = value

"""
    _numerical_model_scalar(expr, parameters)

Evaluate a scalar expression using resolved model-level parameters.
"""
function _numerical_model_scalar(expr::ScalarExpr, parameters)
    if expr.head == :const
        return expr.args[1]
    elseif expr.head == :param
        key = _path_key(_param_parts(expr.args[1]))
        haskey(parameters, key) || error("Model parameter $key is not defined.")
        return parameters[key]
    elseif expr.head == :+
        return reduce(+, (_numerical_model_scalar(arg, parameters) for arg in expr.args))
    elseif expr.head == :*
        return reduce(*, (_numerical_model_scalar(arg, parameters) for arg in expr.args))
    elseif expr.head == :/
        return _numerical_model_scalar(expr.args[1], parameters) /
               _numerical_model_scalar(expr.args[2], parameters)
    elseif expr.head == :^
        return _numerical_model_scalar(expr.args[1], parameters) ^ expr.args[2]
    elseif expr.head == :neg
        return -_numerical_model_scalar(expr.args[1], parameters)
    elseif expr.head == :cos
        return cos(_numerical_model_scalar(expr.args[1], parameters))
    elseif expr.head == :sin
        return sin(_numerical_model_scalar(expr.args[1], parameters))
    end
    error("Scalar expression head $(expr.head) cannot be evaluated in a model.")
end

"""
    _embed_operator(operator, target, order, dimensions)

Tensor an operator with identities so it acts on `target` in the model space.
"""
function _embed_operator(operator, target, order, dimensions)
    factors = Any[
        id == target ? operator : qeye(dimensions[id])
        for id in order
    ]
    return length(factors) == 1 ? factors[1] : tensor(factors...)
end

"""
    _project_operator(states, operator)

Project `operator` as `states' * operator * states` and return the result as a
`QuantumObject`.
"""
function _project_operator(states, operator)
    data = adjoint(states) * operator.data * states
    return QuantumObject(data)
end
