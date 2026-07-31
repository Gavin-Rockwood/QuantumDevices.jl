"""Options used when dressing product states into interacting model states."""
struct DressingSpec
    steps::Int
    schedule::Symbol
    exponent::Float64
    minimum_overlap::Float64
    diagnostics::Symbol
end

function DressingSpec(;
        steps::Integer = 20,
        schedule::Symbol = :power,
        exponent::Real = 3.0,
        minimum_overlap::Real = 0.5,
        diagnostics::Symbol = :none)
    steps >= 2 || error("DressingSpec steps must be at least two.")
    schedule in (:linear, :power) || error("Unknown dressing schedule $schedule.")
    exponent > 0 || error("DressingSpec exponent must be positive.")
    0 <= minimum_overlap <= 1 || error("Minimum overlap must be between zero and one.")
    diagnostics in (:none, :events, :full) || error("Unknown diagnostics mode $diagnostics.")
    DressingSpec(steps, schedule, Float64(exponent), Float64(minimum_overlap), diagnostics)
end

"""A flat, ordered collection of components and their symbolic interactions."""
struct ModelSpec{S<:NamedTuple,I<:Tuple}
    name::Symbol
    subsystems::S
    interactions::I
    parameters::Dict{ParamPath,DeviceParameter}
    defaults::Dict{Tuple{Vararg{Symbol}},Any}
    dimension::Dimension
    states_to_keep::Union{Nothing,Vector{Tuple{Vararg{Int}}}}
    dressingspec::DressingSpec
end

function ModelSpec(
        name::Symbol,
        components::AbstractVector{<:Component},
        dims = Dict();
        interactions = (),
        pretruncation_dims = Dict(),
        parameters = (),
        defaults = Dict(),
        states_to_keep = nothing,
        dressingspec::DressingSpec = DressingSpec())
    names = Tuple(component.name for component in components)
    length(unique(names)) == length(names) || error("Subsystems contain duplicate component names.")
    all(name -> name in names, keys(dims)) || error("dims contains an unknown component name.")
    all(name -> name in names, keys(pretruncation_dims)) ||
        error("pretruncation_dims contains an unknown component name.")
    all(interaction -> interaction isa OperatorExpr, interactions) ||
        error("Interactions must be operator expressions.")
    model_parameters = DeviceParameter[parameters...]
    provided_defaults = Dict{Tuple{Vararg{Symbol}},Any}(
        _model_default_parts(path) => value for (path, value) in defaults
    )

    prepared = Component[
        haskey(pretruncation_dims, component.name) ?
            _pretruncate(component, Int(pretruncation_dims[component.name]), provided_defaults) : component
        for component in components
    ]
    for component in prepared
        length(component.spec.dimension) == 1 ||
            error("Component $(component.name) must have one effective dimension.")
        isfinite(only(component.spec.dimension)) ||
            error("Component $(component.name) has infinite dimension; provide pretruncation_dims[:$(component.name)].")
    end

    final_dims = Int[
        haskey(dims, component.name) ? Int(dims[component.name]) : Int(only(component.spec.dimension))
        for component in prepared
    ]
    all(>(0), final_dims) || error("Model dimensions must be positive.")
    for (component, dimension) in zip(prepared, final_dims)
        dimension <= only(component.spec.dimension) ||
            error("Model requests dimension $dimension from component $(component.name) with dimension $(only(component.spec.dimension)).")
    end

    systems = NamedTuple{names}(Tuple(prepared))
    registry = Dict{ParamPath,DeviceParameter}()
    for (component_name, component) in pairs(systems)
        for (path, parameter) in QuantumDevices.parameters(component)
            registry[ParamPath(component_name, path.parts...)] = parameter
        end
    end
    for parameter in model_parameters
        registry[parameter.path] = parameter
    end
    for interaction in interactions, path in keys(QuantumDevices.parameters(interaction))
        haskey(registry, path) ||
            error("Interaction parameter $(join(path.parts, "/")) was not passed in parameters.")
    end
    for (path, value) in provided_defaults
        parameter_path = ParamPath(path...)
        haskey(registry, parameter_path) ||
            error("Default path $(join(path, "/")) is not a registered model parameter.")
        _parameter_value(registry[parameter_path], value)
    end
    default_values = Dict{Tuple{Vararg{Symbol}},Any}()
    for (path, parameter) in registry
        value = get(provided_defaults, path.parts, parameter.default)
        _parameter_value(parameter, value)
        default_values[path.parts] = value
    end

    selected = states_to_keep === nothing ? nothing :
        Tuple{Vararg{Int}}[Tuple(Int.(state)) for state in states_to_keep]
    selected === nothing || length(unique(selected)) == length(selected) ||
        error("states_to_keep contains duplicate states.")
    if selected !== nothing
        for state in selected
            length(state) == length(final_dims) || error("State $state has the wrong number of indices.")
            all(0 <= state[i] < final_dims[i] for i in eachindex(final_dims)) ||
                error("State $state is outside the model dimension.")
        end
    end

    ModelSpec(name, systems, Tuple(interactions), registry, default_values,
        Dimension(Tuple(final_dims)), selected, dressingspec)
end

_model_default_parts(path::Union{Tuple,AbstractVector}) =
    Tuple(Symbol.(path))
_model_default_parts(path) = _param_path(path).parts

"""Numerical model built from one `ModelSpec`."""
struct QuantumDeviceModel
    spec::ModelSpec
    hamiltonian
    operators::Dict{Tuple{Vararg{Symbol}},Any}
    states::Dict{Tuple{Vararg{Int}},Any}
    energies::Dict{Tuple{Vararg{Int}},Float64}
    dressing_res::StateTrackingResult
end

"""
    parameters(model)

Return the active non-fixed parameters accepted by `model.hamiltonian`.
The dictionary keys are the `ParamPath`s expected by the `QobjEvo` parameter
argument.
"""
parameters(model::QuantumDeviceModel) = _hamiltonian_parameters(model.spec)

"""
    hamiltonian(model; param=Dict{ParamPath,Any}(), t=0.0)

Evaluate a model Hamiltonian at time `t` using `ParamPath`-keyed overrides.
When `model.hamiltonian` is a `QobjEvo`, its direct call interface accepts the
same dictionary as `model.hamiltonian(param, t)`.
"""
function hamiltonian(
        model::QuantumDeviceModel;
        param = Dict{ParamPath,Any}(),
        t::Real = 0.0)
    active = parameters(model)
    checked = _validated_hamiltonian_param(model.spec, active, param)
    model.hamiltonian isa QuantumObjectEvolution ?
        model.hamiltonian(checked, t) : model.hamiltonian
end

"""Evaluate a symbolic expression in a built model's operator basis."""
function numerical(model::QuantumDeviceModel, expression::OperatorExpr)
    _numerical_expression(
        expression,
        path -> begin
            key = _canonical_reference(path, :operators)
            haskey(model.operators, key) ||
                error("Operator path $(join(key, "/")) is not available in model $(model.spec.name).")
            model.operators[key]
        end,
        parameter -> _model_value(model.spec, _model_parameter_path(model.spec, parameter)),
    )
end

numerical(model::QuantumDeviceModel, expression::ScalarExpr) =
    _numerical_scalar(expression,
        parameter -> _model_value(model.spec, _model_parameter_path(model.spec, parameter)))

"""
    model(component::Component; dim=nothing, params=nothing)

Build a temporary one-component model. `params` may be any mapping-like object
whose keys are local component parameter paths.
"""
function model(component::Component; dim = nothing, params = nothing)
    length(component.spec.dimension) == 1 ||
        error("Component $(component.name) must have one effective dimension.")
    maximum_dimension = only(component.spec.dimension)
    if dim === nothing
        isfinite(maximum_dimension) ||
            error("Component $(component.name) has infinite dimension; provide dim.")
        dimension = Int(maximum_dimension)
    else
        dimension = Int(only(Dimension(dim)))
    end

    dims = Dict(component.name => dimension)
    pretruncation_dims = isfinite(maximum_dimension) ? Dict{Symbol,Int}() : copy(dims)
    defaults = _component_model_defaults(component, params)
    spec = ModelSpec(
        component.name,
        Component[component],
        dims;
        pretruncation_dims,
        defaults,
    )
    model(spec)
end

function _component_model_defaults(component, params)
    Dict{Tuple{Vararg{Symbol}},Any}(
        (component.name, path.parts...) => value
        for (path, value) in _component_parameter_overrides(component, params)
    )
end

function model(spec::ModelSpec)
    names = keys(spec.subsystems)
    dims = spec.dimension.size
    components = values(spec.subsystems)
    local_hamiltonians = Any[]
    local_operators = Dict{Tuple{Vararg{Symbol}},Any}()

    for (name, component, dimension) in zip(names, components, dims)
        values = Dict(
            _local_key(path.parts[2:end]) => _at(_model_value(spec, path), 0)
            for path in keys(spec.parameters)
            if length(path.parts) >= 2 && path.parts[1] == name
        )
        push!(local_hamiltonians, numerical(component; dimension = (dimension,), params = values))
        for operator in operators(component)
            local_operator = numerical(component, operator; dimension = (dimension,), params = values)
            local_operators[(name, _operator_parts(operator)...)] = local_operator
        end
    end

    embedded_operators = Dict(path => _embed(operator, path[1], names, dims)
        for (path, operator) in local_operators)
    uncoupled = sum(_embed(hamiltonian, name, names, dims)
        for (hamiltonian, name) in zip(local_hamiltonians, names))
    interaction = isempty(spec.interactions) ? 0 * uncoupled : sum(
        _numerical_expression(expr,
            path -> embedded_operators[_canonical_reference(path, :operators)],
            parameter -> _at(_model_value(spec, _model_parameter_path(spec, parameter)), 0))
        for expr in spec.interactions
    )

    labels = if spec.states_to_keep === nothing
        vec([Tuple(state) for state in Iterators.product(
            (0:(dimension - 1) for dimension in dims)...,
        )])
    else
        spec.states_to_keep
    end
    local_states = [eigenstates(hamiltonian).vectors for hamiltonian in local_hamiltonians]
    references = hcat((foldl(kron,
        (local_states[i][:, state[i] + 1] for i in eachindex(dims)); init = ComplexF64[1])
        for state in labels)...)

    dressing = spec.dressingspec
    dressing_res = get_dressed_states(uncoupled, interaction, references;
        steps = dressing.steps,
        schedule = dressing.schedule,
        exponent = dressing.exponent,
        minimum_overlap = dressing.minimum_overlap,
        diagnostics = dressing.diagnostics)
    dressed_states = hcat((_state_vector(state) for state in dressing_res.states[:, end])...)
    dressed_energies = Float64.(real.(dressing_res.other_sorts[:energy][:, end]))
    states = Dict(label => QuantumObject(dressed_states[:, i]) for (i, label) in enumerate(labels))
    energies = Dict(label => dressed_energies[i] for (i, label) in enumerate(labels))

    selected = spec.states_to_keep !== nothing
    default_hamiltonian = selected ?
        QuantumObject(Diagonal(ComplexF64.(dressed_energies))) :
        uncoupled + interaction
    final_operators = selected ?
        Dict(path => _project(dressed_states, operator) for (path, operator) in embedded_operators) : embedded_operators
    stored_hamiltonian =
        _compiled_model_hamiltonian(spec, final_operators, default_hamiltonian)
    QuantumDeviceModel(
        spec,
        stored_hamiltonian,
        final_operators,
        states,
        energies,
        dressing_res,
    )
end

function _hamiltonian_parameters(spec::ModelSpec)
    registry = Dict{ParamPath,DeviceParameter}()
    for (name, component) in pairs(spec.subsystems)
        for (path, _) in parameters(component)
            model_path = ParamPath(name, path.parts...)
            parameter = spec.parameters[model_path]
            parameter.fixed || (registry[model_path] = parameter)
        end
    end
    for interaction in spec.interactions
        for parameter in values(parameters(interaction))
            path = _model_parameter_path(spec, parameter)
            registered = spec.parameters[path]
            registered.fixed || (registry[path] = registered)
        end
    end
    registry
end

function _validated_hamiltonian_param(spec, active, param)
    param === nothing && return Dict{ParamPath,Any}()
    param isa AbstractDict ||
        error("Hamiltonian param must be an AbstractDict with ParamPath keys.")
    for path in keys(param)
        path isa ParamPath ||
            error("Hamiltonian param keys must be ParamPath objects; got $(typeof(path)).")
        if !haskey(active, path)
            if haskey(spec.parameters, path) && spec.parameters[path].fixed
                error("Hamiltonian parameter $(join(path.parts, "/")) is fixed.")
            end
            error("Hamiltonian parameter $(join(path.parts, "/")) is not an active non-fixed parameter.")
        end
    end
    param
end

function _compiled_model_hamiltonian(spec, operators, default_hamiltonian)
    _compile_model_hamiltonian(spec, operators, default_hamiltonian, nothing)
end

function _bound_model_hamiltonian(model::QuantumDeviceModel, controls)
    _compile_model_hamiltonian(
        model.spec,
        model.operators,
        hamiltonian(model),
        controls,
    )
end

function _compile_model_hamiltonian(spec, operators, default_hamiltonian, controls)
    active = _hamiltonian_parameters(spec)
    isempty(active) && return default_hamiltonian

    total = nothing
    for (name, component) in pairs(spec.subsystems)
        value = _numerical_expression(
            hamiltonian(component),
            path -> _model_operator(operators, path, name),
            parameter -> begin
                path = ParamPath(name, parameter.path.parts...)
                _model_hamiltonian_coefficient(spec, active, path, controls)
            end,
        )
        total = total === nothing ? value : total + value
    end
    for interaction in spec.interactions
        value = _numerical_expression(
            interaction,
            path -> _model_operator(operators, path),
            parameter -> _model_hamiltonian_coefficient(
                spec,
                active,
                _model_parameter_path(spec, parameter),
                controls,
            ),
        )
        total = total === nothing ? value : total + value
    end

    baseline = total isa QuantumObjectEvolution ?
        total(Dict{ParamPath,Any}(), 0.0) : total
    total + default_hamiltonian - baseline
end

function _model_hamiltonian_coefficient(spec, active, path, controls = nothing)
    parameter = spec.parameters[path]
    default = _model_value(spec, path)
    parameter.fixed && return default
    if controls !== nothing
        value = get(controls, path, default)
        return _parameter_value(parameter, value)
    end
    _ParameterizedCoefficient((param, time) -> begin
        checked = _validated_hamiltonian_param(spec, active, param)
        value = haskey(checked, path) ?
            _parameter_value(parameter, checked[path]) :
            default
        _at(value, time)
    end)
end

function _model_operator(operators, path, component = nothing)
    rooted = length(path) >= 1 && path[1] == :components
    relative = _canonical_reference(path, :operators)
    key = rooted || component === nothing ? relative : (component, relative...)
    haskey(operators, key) ||
        error("Operator path $(join(key, "/")) is not available in model Hamiltonian.")
    operators[key]
end

function _model_parameter_path(spec::ModelSpec, parameter::DeviceParameter)
    matches = ParamPath[path for (path, candidate) in spec.parameters if candidate === parameter]
    length(matches) == 1 && return only(matches)
    haskey(spec.parameters, parameter.path) && return parameter.path
    isempty(matches) &&
        error("Parameter $(join(parameter.path.parts, "/")) is not registered in model $(spec.name).")
    error("Parameter $(join(parameter.path.parts, "/")) has multiple paths in model $(spec.name).")
end

function _model_value(spec::ModelSpec, path::ParamPath)
    haskey(spec.parameters, path) ||
        error("Parameter $(join(path.parts, "/")) is not registered in model $(spec.name).")
    parameter = spec.parameters[path]
    haskey(spec.defaults, path.parts) ||
        error("ModelSpec $(spec.name) does not snapshot parameter $(join(path.parts, "/")).")
    _parameter_value(parameter, spec.defaults[path.parts])
end

function _pretruncate(component::Component, dimension::Int, defaults)
    dimension > 0 || error("Pretruncation dimensions must be positive.")
    values = Dict(
        _local_key(path.parts) => _at(_parameter_value(
            parameter,
            get(defaults, (component.name, path.parts...), parameter.default),
        ), 0)
        for (path, parameter) in parameters(component)
    )
    local_hamiltonian = numerical(component; dimension = (dimension,), params = values)
    size(local_hamiltonian, 1) == dimension ||
        error("Component $(component.name) produced dimension $(size(local_hamiltonian, 1)), expected $dimension.")
    eigensystem = eigenstates(local_hamiltonian)
    projection = eigensystem.vectors
    projected_operators = Dict(operator => _projected_operator(
        _project(projection, numerical(component, operator; dimension = (dimension,), params = values)),
    ) for operator in operators(component))
    spec = GenericSpec(real.(eigensystem.values);
        operators = projected_operators,
        hamiltonian = hamiltonian(component),
        source = component.spec)
    Component(spec, component.name)
end

function component(model::QuantumDeviceModel; name = model.spec.name, dimension = nothing)
    eigensystem = eigenstates(hamiltonian(model))
    keep = dimension === nothing ? length(eigensystem.values) : Int(dimension)
    1 <= keep <= length(eigensystem.values) || error("Invalid effective component dimension $keep.")
    projection = eigensystem.vectors[:, 1:keep]
    projected_operators = Dict(path => _projected_operator(_project(projection, operator))
        for (path, operator) in model.operators)
    Component(GenericSpec(real.(eigensystem.values[1:keep]);
        operators = projected_operators, source = model.spec), Symbol(name))
end

_operator_parts(operator::Symbol) = (operator,)
_operator_parts(operator::Tuple) = operator
_local_key(path::Tuple) = length(path) == 1 ? only(path) : path

function _embed(operator, target, names, dims)
    factors = Any[name == target ? operator : qeye(dimension) for (name, dimension) in zip(names, dims)]
    length(factors) == 1 ? only(factors) : tensor(factors...)
end

_project(states, operator) = QuantumObject(adjoint(states) * operator.data * states)

_projected_operator(operator) = (; dimension, kwargs...) -> begin
    keep = only(dimension)
    keep <= size(operator, 1) || error("Projected operator has dimension $(size(operator, 1)), received $keep.")
    QuantumObject(operator.data[1:keep, 1:keep])
end
