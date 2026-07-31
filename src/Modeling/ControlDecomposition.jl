struct ModelControlTerm{F}
    operator::Matrix{ComplexF64}
    coefficient::F
    linear_weights::Union{Nothing,Vector{Float64}}
end

struct ModelControlDecomposition
    controls::Vector{ParamPath}
    drift::Matrix{ComplexF64}
    terms::Vector{ModelControlTerm}
end

function _model_control_decomposition(
    model::QuantumDeviceModel,
    control_paths::Vector{ParamPath},
)
    active = parameters(model)
    isempty(control_paths) &&
        error("Model control decomposition requires at least one control.")
    length(unique(control_paths)) == length(control_paths) ||
        error("Model control paths must be unique.")
    for path in control_paths
        haskey(active, path) ||
            error("Model control $(join(path.parts, "/")) is not active.")
    end
    indices = Dict(path => index for (index, path) in enumerate(control_paths))
    default_controls = Float64[]
    for path in control_paths
        value = _at(_model_value(model.spec, path), 0.0)
        value isa Real || error("Model control $(join(path.parts, "/")) must be real.")
        push!(default_controls, Float64(value))
    end

    raw_terms = Tuple{Any,Any}[]
    for (name, component) in pairs(model.spec.subsystems)
        append!(
            raw_terms,
            _operator_terms(
                hamiltonian(component),
                path -> _model_operator(model.operators, path, name),
                parameter -> begin
                    path = ParamPath(name, parameter.path.parts...)
                    _decomposition_coefficient(model.spec, indices, path)
                end,
            ),
        )
    end
    for interaction in model.spec.interactions
        append!(
            raw_terms,
            _operator_terms(
                interaction,
                path -> _model_operator(model.operators, path),
                parameter -> _decomposition_coefficient(
                    model.spec,
                    indices,
                    _model_parameter_path(model.spec, parameter),
                ),
            ),
        )
    end

    dimension = size(model.hamiltonian, 1)
    raw_default = zeros(ComplexF64, dimension, dimension)
    drift = zeros(ComplexF64, dimension, dimension)
    terms = ModelControlTerm[]
    for (coefficient, operator) in raw_terms
        matrix = ComplexF64.(Matrix(operator.data))
        value = coefficient isa Function ? coefficient(default_controls) : coefficient
        raw_default .+= value .* matrix
        if coefficient isa Function
            _append_hermitian_control_terms!(
                terms,
                matrix,
                coefficient,
                length(control_paths),
            )
        else
            drift .+= _hermitian_part(coefficient .* matrix)
        end
    end
    drift .+= ComplexF64.(Matrix(hamiltonian(model).data)) .- _hermitian_part(raw_default)

    isapprox(drift, adjoint(drift); atol = 1e-10, rtol = 1e-10) ||
        error("Model control decomposition produced a non-Hermitian drift.")
    decomposition = ModelControlDecomposition(control_paths, drift, terms)
    _validate_control_decomposition(model, decomposition, default_controls)
    decomposition
end

function _decomposition_coefficient(spec, indices, path)
    haskey(indices, path) && return controls -> controls[indices[path]]
    _model_value(spec, path)
end

function _append_hermitian_control_terms!(terms, operator, coefficient, count)
    symmetric = _hermitian_part(operator)
    antisymmetric = (operator - adjoint(operator)) / (2im)
    if norm(symmetric) > 1e-12
        real_coefficient = controls -> real(coefficient(controls))
        push!(
            terms,
            ModelControlTerm(
                symmetric,
                real_coefficient,
                _linear_weights(real_coefficient, count),
            ),
        )
    end
    if norm(antisymmetric) > 1e-12
        imaginary_coefficient = controls -> -imag(coefficient(controls))
        push!(
            terms,
            ModelControlTerm(
                ComplexF64.(antisymmetric),
                imaginary_coefficient,
                _linear_weights(imaginary_coefficient, count),
            ),
        )
    end
end

_hermitian_part(matrix) = (matrix + adjoint(matrix)) / 2

function _linear_weights(coefficient, count)
    origin = zeros(Float64, count)
    offset = coefficient(origin)
    isapprox(offset, 0; atol = 1e-12, rtol = 1e-12) || return nothing
    weights = Float64[]
    for index = 1:count
        basis = zeros(Float64, count)
        basis[index] = 1.0
        push!(weights, Float64(coefficient(basis)))
    end
    probes =
        ([0.173 * index for index = 1:count], [-0.271 + 0.119 * index for index = 1:count])
    all(
        probe ->
            isapprox(coefficient(probe), dot(weights, probe); atol = 1e-10, rtol = 1e-10),
        probes,
    ) ? weights : nothing
end

function _control_decomposition_matrix(decomposition, controls)
    value = copy(decomposition.drift)
    for term in decomposition.terms
        value .+= term.coefficient(controls) .* term.operator
    end
    value
end

function _validate_control_decomposition(model, decomposition, defaults)
    probes = (defaults, defaults .+ [0.013 * index for index in eachindex(defaults)])
    active = parameters(model)
    for probe in probes
        valid = true
        for (path, value) in zip(decomposition.controls, probe)
            try
                _parameter_value(active[path], value)
            catch
                valid = false
                break
            end
        end
        valid || continue
        param = Dict(path => value for (path, value) in zip(decomposition.controls, probe))
        expected = hamiltonian(model; param, t = 0.0)
        actual = _control_decomposition_matrix(decomposition, probe)
        isapprox(actual, Matrix(expected.data); atol = 1e-8, rtol = 1e-8) ||
            error("Model control decomposition failed reconstruction validation.")
    end
    decomposition
end

"""A Piccolo pulse-optimization request over selected `ModelSpec` controls."""
struct PiccoloOptimizationProblem{O<:AbstractGateObjective}
    name::Symbol
    modelspec::ModelSpec
    controls::Vector{ParamPath}
    bounds::Vector{Tuple{Float64,Float64}}
    duration::Float64
    knots::Int
    initial_controls::Matrix{Float64}
    control_modulations::Vector{Any}
    objective::O
    template_kwargs::Dict{Symbol,Any}
    solve_kwargs::Dict{Symbol,Any}
end

function PiccoloOptimizationProblem(
    name::Symbol,
    modelspec::ModelSpec,
    bounds::AbstractDict,
    objective::AbstractGateObjective;
    duration::Real,
    knots::Integer = 100,
    initial_controls = nothing,
    control_modulations = Dict{ParamPath,Any}(),
    template_kwargs = Dict{Symbol,Any}(),
    solve_kwargs = Dict{Symbol,Any}(),
)
    duration > 0 || error("Piccolo gate duration must be positive.")
    knots >= 3 || error("Piccolo optimization requires at least three knots.")
    all(path -> path isa ParamPath, keys(bounds)) ||
        error("Piccolo control bounds must use ParamPath keys.")
    control_modulations isa AbstractDict ||
        error("control_modulations must be an AbstractDict.")
    all(path -> path isa ParamPath, keys(control_modulations)) ||
        error("Piccolo control modulation keys must be ParamPaths.")
    all(path -> haskey(bounds, path), keys(control_modulations)) ||
        error("Piccolo control modulations must refer to selected control paths.")
    controls = sort!(ParamPath[keys(bounds)...]; by = path -> path.parts)
    active = _hamiltonian_parameters(modelspec)
    prepared_bounds = Tuple{Float64,Float64}[]
    defaults = Float64[]
    prepared_modulations = Any[]
    for path in controls
        haskey(active, path) ||
            error("Piccolo control $(join(path.parts, "/")) is not active.")
        bound = bounds[path]
        bound isa Tuple && length(bound) == 2 ||
            error("Piccolo bounds for $(join(path.parts, "/")) must be a pair.")
        lower, upper = Float64.(bound)
        lower < upper ||
            error("Piccolo bounds for $(join(path.parts, "/")) must be ordered.")
        default = _at(_model_value(modelspec, path), 0.0)
        default isa Real || error("Piccolo controls must have real defaults.")
        lower <= default <= upper ||
            error("Default for $(join(path.parts, "/")) lies outside its bounds.")
        _parameter_value(active[path], lower)
        _parameter_value(active[path], upper)
        push!(prepared_bounds, (lower, upper))
        push!(defaults, Float64(default))
        modulation = get(control_modulations, path, nothing)
        if modulation !== nothing
            applicable(modulation, 0.0) || error(
                "Piccolo modulation for $(join(path.parts, "/")) must accept time.",
            )
            all(
                time -> begin
                    value = modulation(time)
                    value isa Real && isfinite(value)
                end,
                (0.0, Float64(duration) / 2, Float64(duration)),
            ) || error(
                "Piccolo modulation for $(join(path.parts, "/")) must return finite real values.",
            )
            iszero(default) || error(
                "Carrier-modulated Piccolo control $(join(path.parts, "/")) must have a zero model default.",
            )
        end
        push!(prepared_modulations, modulation)
    end
    isempty(controls) && error("Piccolo optimization requires at least one control.")

    prepared_initial = if initial_controls === nothing
        repeat(defaults, 1, Int(knots))
    else
        Matrix{Float64}(initial_controls)
    end
    size(prepared_initial) == (length(controls), Int(knots)) ||
        error("initial_controls must have size " * "($(length(controls)), $(Int(knots))).")
    for (index, bound) in enumerate(prepared_bounds)
        all(value -> bound[1] <= value <= bound[2], prepared_initial[index, :]) ||
            error("Initial Piccolo controls violate bounds for $(controls[index]).")
        prepared_initial[index, 1] = defaults[index]
        prepared_initial[index, end] = defaults[index]
    end
    PiccoloOptimizationProblem(
        name,
        modelspec,
        controls,
        prepared_bounds,
        Float64(duration),
        Int(knots),
        prepared_initial,
        prepared_modulations,
        objective,
        Dict{Symbol,Any}(template_kwargs),
        Dict{Symbol,Any}(solve_kwargs),
    )
end
