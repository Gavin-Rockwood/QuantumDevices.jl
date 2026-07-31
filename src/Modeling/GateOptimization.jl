"""One named, bounded scalar used by an envelope gate builder."""
struct GateVariable
    name::Symbol
    initial::Float64
    lower::Float64
    upper::Float64
end

function GateVariable(name::Symbol; initial::Real, lower::Real, upper::Real)
    lower < upper || error("GateVariable $name must have lower < upper.")
    lower <= initial <= upper ||
        error("GateVariable $name initial value must lie within its bounds.")
    GateVariable(name, Float64(initial), Float64(lower), Float64(upper))
end

abstract type AbstractGateObjective end

"""A phase-insensitive logical-unitary target with optional leakage penalty."""
struct UnitaryObjective{T<:AbstractMatrix,S<:Tuple,F} <: AbstractGateObjective
    target::T
    states::S
    leakage_weight::Float64
    frame_frequencies::F
end

function UnitaryObjective(
    target;
    states,
    leakage_weight::Real = 1.0,
    frame_frequencies = nothing,
)
    labels = Tuple(states)
    isempty(labels) && error("UnitaryObjective states must not be empty.")
    leakage_weight >= 0 || error("UnitaryObjective leakage_weight must be nonnegative.")
    matrix = _objective_matrix(target)
    size(matrix) == (length(labels), length(labels)) || error(
        "UnitaryObjective target has size $(size(matrix)); expected " *
        "($(length(labels)), $(length(labels))).",
    )
    frequencies = _prepare_frame_frequencies(frame_frequencies)
    UnitaryObjective(matrix, labels, Float64(leakage_weight), frequencies)
end

function _prepare_frame_frequencies(frequencies)
    frequencies === nothing && return nothing
    frequencies isa NamedTuple ||
        error("frame_frequencies must be a NamedTuple keyed by subsystem name.")
    isempty(frequencies) && error("frame_frequencies must not be empty.")
    for (name, frequency) in pairs(frequencies)
        frequency isa Real && isfinite(frequency) ||
            error("Frame frequency $name must be a finite real value.")
    end
    NamedTuple{propertynames(frequencies)}(Tuple(Float64.(values(frequencies))))
end

"""An initial-to-target pure-state objective, optionally within a labeled subspace."""
struct StateTransferObjective{I,T,S} <: AbstractGateObjective
    initial::I
    target::T
    states::S
end

function StateTransferObjective(initial, target; states = nothing)
    labels = states === nothing ? nothing : Tuple(states)
    labels === nothing ||
        !isempty(labels) ||
        error("StateTransferObjective states must not be empty.")
    StateTransferObjective(initial, target, labels)
end

"""A bounded classical optimization over a function that constructs `GateSpec`s."""
struct EnvelopeOptimizationProblem{B,O<:AbstractGateObjective}
    builder::B
    variables::Vector{GateVariable}
    objective::O
    solver_kwargs::Dict{Symbol,Any}
end

function EnvelopeOptimizationProblem(
    builder,
    variables,
    objective::AbstractGateObjective;
    solver_kwargs = Dict{Symbol,Any}(),
)
    prepared = GateVariable[variables...]
    isempty(prepared) &&
        error("EnvelopeOptimizationProblem requires at least one GateVariable.")
    names = getfield.(prepared, :name)
    length(unique(names)) == length(names) ||
        error("EnvelopeOptimizationProblem variable names must be unique.")
    solver_kwargs isa AbstractDict || error("solver_kwargs must be an AbstractDict.")
    EnvelopeOptimizationProblem(
        builder,
        prepared,
        objective,
        Dict{Symbol,Any}(solver_kwargs),
    )
end

"""The backend-independent result of classical or optimal-control gate synthesis."""
struct GateOptimizationResult{G,R}
    gate::G
    minimizer::Dict{Symbol,Float64}
    loss::Float64
    fidelity::Float64
    leakage::Float64
    converged::Bool
    evaluations::Int
    settings::Dict{Symbol,Any}
    failures::Vector{NamedTuple}
    backend_result::R
end

function Base.show(io::IO, result::GateOptimizationResult)
    values = join(
        (
            "$name=$(result.minimizer[name])" for
            name in sort!(collect(keys(result.minimizer)))
        ),
        ", ",
    )
    print(
        io,
        "GateOptimizationResult(gate=",
        result.gate === nothing ? "nothing" : result.gate.name,
        isempty(values) ? "" : ", values=($values)",
        ", loss=",
        result.loss,
        ", fidelity=",
        result.fidelity,
        ", leakage=",
        result.leakage,
        ", converged=",
        result.converged,
        ", evaluations=",
        result.evaluations,
        ")",
    )
end

"""
    optimize_gate(problem::EnvelopeOptimizationProblem; kwargs...)

Evaluate deterministic coordinate seeds, then refine the best valid candidate
with bounded derivative-free optimization.
"""
function optimize_gate(
    problem::EnvelopeOptimizationProblem;
    method = Optim.Fminbox(Optim.NelderMead()),
    options = nothing,
    iterations = nothing,
    coarse_fractions = (0.25, 0.5, 0.75),
    show_trace = nothing,
    store_trace = nothing,
)
    prepared_options, option_settings =
        _optim_options(options; iterations, show_trace, store_trace)
    all(fraction -> 0 <= fraction <= 1, coarse_fractions) ||
        error("coarse_fractions must lie between zero and one.")

    variables = problem.variables
    lower = getfield.(variables, :lower)
    upper = getfield.(variables, :upper)
    initial = getfield.(variables, :initial)
    cache = Dict{Tuple{Vararg{Float64}},Any}()
    failures = NamedTuple[]

    function evaluate(values)
        key = Tuple(clamp.(Float64.(values), lower, upper))
        haskey(cache, key) && return cache[key]
        named = NamedTuple{Tuple(getfield.(variables, :name))}(key)
        evaluated = try
            gate = problem.builder(named)
            gate isa GateSpec ||
                error("EnvelopeOptimizationProblem builder must return a GateSpec.")
            metrics = _evaluate_gate_objective(
                gate,
                problem.objective;
                solver_kwargs = problem.solver_kwargs,
            )
            all(isfinite, (metrics.loss, metrics.fidelity, metrics.leakage)) ||
                error("Gate objective returned a nonfinite metric.")
            (; valid = true, gate, metrics..., error = nothing)
        catch exception
            failure = (;
                variables = Dict(
                    variable.name => key[index] for
                    (index, variable) in enumerate(variables)
                ),
                exception,
            )
            push!(failures, failure)
            (;
                valid = false,
                gate = nothing,
                loss = Inf,
                fidelity = 0.0,
                leakage = Inf,
                error = exception,
            )
        end
        cache[key] = evaluated
        evaluated
    end

    seeds = Vector{Vector{Float64}}([copy(initial)])
    for index in eachindex(variables), fraction in coarse_fractions
        candidate = copy(initial)
        candidate[index] = lower[index] + Float64(fraction) * (upper[index] - lower[index])
        any(seed -> seed == candidate, seeds) || push!(seeds, candidate)
    end
    foreach(evaluate, seeds)
    valid_seeds = filter(seed -> evaluate(seed).valid, seeds)

    backend_result = nothing
    backend_converged = false
    if !isempty(valid_seeds)
        seed = copy(argmin(seed -> evaluate(seed).loss, valid_seeds))
        objective_function(values) = begin
            value = evaluate(values).loss
            isfinite(value) ? value : 1e100
        end
        try
            backend_result = Optim.optimize(
                objective_function,
                lower,
                upper,
                seed,
                method,
                prepared_options,
            )
            backend_converged = Optim.converged(backend_result)
            evaluate(Optim.minimizer(backend_result))
        catch exception
            backend_result = exception
            push!(failures, (; variables = Dict{Symbol,Float64}(), exception))
        end
    end

    valid = [(values, evaluation) for (values, evaluation) in cache if evaluation.valid]
    settings = Dict{Symbol,Any}(
        :backend => :Optim,
        :method => method,
        :options => prepared_options,
        :option_settings => option_settings,
        :coarse_fractions => Tuple(Float64.(coarse_fractions)),
        :solver_kwargs => copy(problem.solver_kwargs),
    )
    if isempty(valid)
        return GateOptimizationResult(
            nothing,
            Dict(variable.name => variable.initial for variable in variables),
            Inf,
            0.0,
            Inf,
            false,
            length(cache),
            settings,
            failures,
            backend_result,
        )
    end

    best_values, best = argmin(entry -> entry[2].loss, valid)
    minimizer = Dict(
        variable.name => best_values[index] for (index, variable) in enumerate(variables)
    )
    GateOptimizationResult(
        best.gate,
        minimizer,
        best.loss,
        best.fidelity,
        best.leakage,
        backend_converged,
        length(cache),
        settings,
        failures,
        backend_result,
    )
end

function _optim_options(
    options;
    iterations = nothing,
    show_trace = nothing,
    store_trace = nothing,
)
    convenience = (; iterations, show_trace, store_trace)
    if options !== nothing
        options isa Optim.Options || error("options must be an Optim.Options object.")
        all(isnothing, values(convenience)) || error(
            "Pass either options=Optim.Options(...) or direct optimizer options, not both.",
        )
        return options, Dict{Symbol,Any}(:source => :options)
    end
    prepared_iterations = iterations === nothing ? 200 : Int(iterations)
    prepared_iterations > 0 || error("iterations must be positive.")
    prepared_show_trace = show_trace === nothing ? false : Bool(show_trace)
    prepared_store_trace = store_trace === nothing ? false : Bool(store_trace)
    prepared = Optim.Options(
        iterations = prepared_iterations,
        show_trace = prepared_show_trace,
        store_trace = prepared_store_trace,
    )
    prepared,
    Dict{Symbol,Any}(
        :iterations => prepared_iterations,
        :show_trace => prepared_show_trace,
        :store_trace => prepared_store_trace,
    )
end

# Implemented by QuantumDevicesPiccoloExt when Piccolo is loaded.
function _optimize_piccolo end

"""
    optimize(gate, target, parameters_to_optimize; kwargs...)

Optimize a parameterized `GateSpec` without mutating it. Inline selector bounds
may be supplied as one pair or a tuple of pairs; selectors may instead be paired
with a `bounds` dictionary or named tuple.
"""
function optimize(
    gate::GateSpec,
    target,
    parameters_to_optimize;
    backend::Symbol = :optim,
    bounds = nothing,
    states = nothing,
    frame_frequencies = nothing,
    leakage_weight::Real = 1.0,
    solver_kwargs = Dict{Symbol,Any}(),
    method = nothing,
    options = nothing,
    iterations = nothing,
    show_trace = nothing,
    store_trace = nothing,
    coarse_fractions = (0.25, 0.5, 0.75),
    kwargs...,
)
    selectors, prepared_bounds = _optimization_selectors(parameters_to_optimize, bounds)
    objective = _optimization_objective(
        gate,
        target,
        states,
        leakage_weight,
        frame_frequencies,
    )
    if backend === :optim
        if !isempty(kwargs)
            suffix = length(kwargs) == 1 ? "" : "s"
            error("Unsupported Optim keyword$suffix: " * join(string.(keys(kwargs)), ", "))
        end
        all(selector -> selector isa Symbol, selectors) ||
            error("The Optim backend accepts only named GateSpec parameters and :duration.")
        gate.control_recipes === nothing && error(
            "The Optim backend requires a GateSpec with named parameters and a control recipe.",
        )
        variables = GateVariable[]
        for selector in selectors
            selector === :duration ||
                hasproperty(gate.parameters, selector) ||
                error("GateSpec $(gate.name) has no named parameter $selector.")
            initial = selector === :duration ? gate.duration : gate.parameters[selector]
            lower, upper = prepared_bounds[selector]
            push!(variables, GateVariable(selector; initial, lower, upper))
        end
        builder = values -> with_parameters(gate; pairs(values)...)
        problem = EnvelopeOptimizationProblem(builder, variables, objective; solver_kwargs)
        prepared_method = method === nothing ? Optim.Fminbox(Optim.NelderMead()) : method
        return optimize_gate(
            problem;
            method = prepared_method,
            options,
            iterations,
            show_trace,
            store_trace,
            coarse_fractions,
        )
    elseif backend === :piccolo
        method === nothing || error(
            "method is an Optim-only keyword and cannot be used with backend=:piccolo.",
        )
        all(isnothing, (options, iterations, show_trace, store_trace)) || error(
            "Optim options cannot be used with backend=:piccolo; use solve_kwargs instead.",
        )
        Base.get_extension(@__MODULE__, :QuantumDevicesPiccoloExt) === nothing && error(
            "Piccolo optimization is unavailable. Install Piccolo.jl and load it with `using Piccolo` before calling optimize(...; backend=:piccolo).",
        )
        return _optimize_piccolo(
            gate,
            objective,
            selectors,
            prepared_bounds;
            solver_kwargs,
            kwargs...,
        )
    end
    error("Unknown gate optimization backend :$backend; expected :optim or :piccolo.")
end

function _optimization_selectors(input, supplied_bounds)
    entries =
        input isa Pair ? (input,) :
        input isa Union{Tuple,AbstractVector} ? Tuple(input) : (input,)
    isempty(entries) && error("At least one gate parameter must be selected.")
    inline = all(entry -> entry isa Pair, entries)
    any(entry -> entry isa Pair, entries) &&
        !inline &&
        error("Optimization selectors cannot mix bounded pairs and bare selectors.")
    inline &&
        supplied_bounds !== nothing &&
        error("Bounds were supplied both inline and through bounds=.")
    selectors = Any[inline ? first(entry) : entry for entry in entries]
    length(unique(selectors)) == length(selectors) ||
        error("Optimization selectors must be unique.")
    all(selector -> selector isa Symbol || selector isa ParamPath, selectors) ||
        error("Optimization selectors must be Symbols or ParamPaths.")
    bound_values =
        inline ? Any[last(entry) for entry in entries] :
        begin
            supplied_bounds === nothing &&
                error("Bare optimization selectors require bounds=.")
            [_optimization_bound(supplied_bounds, selector) for selector in selectors]
        end
    prepared = Dict{Any,Tuple{Float64,Float64}}()
    for (selector, value) in zip(selectors, bound_values)
        (value isa Tuple || value isa Pair) && length(value) == 2 ||
            error("Bounds for $selector must contain exactly two values.")
        lower, upper = Float64.(Tuple(value))
        isfinite(lower) && isfinite(upper) && lower < upper ||
            error("Bounds for $selector must be finite and ordered.")
        prepared[selector] = (lower, upper)
    end
    selectors, prepared
end

function _optimization_bound(bounds, selector)
    if bounds isa NamedTuple
        selector isa Symbol && hasproperty(bounds, selector) ||
            error("bounds does not contain selector $selector.")
        return bounds[selector]
    elseif bounds isa AbstractDict
        haskey(bounds, selector) || error("bounds does not contain selector $selector.")
        return bounds[selector]
    end
    error("bounds must be a NamedTuple or AbstractDict.")
end

function _optimization_objective(
    gate,
    target::AbstractGateObjective,
    states,
    leakage_weight,
    frame_frequencies,
)
    states === nothing || error("states is already defined by the supplied gate objective.")
    leakage_weight == 1 ||
        error("leakage_weight is already defined by the supplied gate objective.")
    frame_frequencies === nothing ||
        error("frame_frequencies is already defined by the supplied gate objective.")
    target
end

function _optimization_objective(gate, target, states, leakage_weight, frame_frequencies)
    matrix = _objective_matrix(target)
    labels = states === nothing ? _infer_gate_states(gate, size(matrix, 1)) : Tuple(states)
    UnitaryObjective(matrix; states = labels, leakage_weight, frame_frequencies)
end

function _infer_gate_states(gate, target_dimension)
    spec = gate.modelspec
    if spec.states_to_keep !== nothing && length(spec.states_to_keep) == target_dimension
        return Tuple(spec.states_to_keep)
    end
    built = model(spec)
    if length(spec.subsystems) == 1 && size(built.hamiltonian, 1) == target_dimension
        return Tuple(sort!(collect(keys(built.states)); by = identity))
    end
    error(
        "Cannot infer an unambiguous logical-state order for target dimension " *
        "$target_dimension; pass states=[...] explicitly.",
    )
end

function _evaluate_gate_objective(
    gate::GateSpec,
    objective::AbstractGateObjective;
    solver_kwargs = Dict{Symbol,Any}(),
)
    built = model(gate.modelspec)
    propagator = _gate_propagator(built, gate; solver_kwargs)
    _objective_metrics(built, propagator, objective, gate.duration)
end

function _gate_propagator(model, gate; solver_kwargs)
    gate_hamiltonian = numerical(model, gate)
    kwargs = Dict{Symbol,Any}(:progress_bar => Val(false))
    merge!(kwargs, solver_kwargs)
    initial = qeye_like(gate_hamiltonian)
    initial = initial isa QuantumObjectEvolution ? initial(nothing, 0.0) : initial
    solution = sesolve(2π * gate_hamiltonian, initial, [0.0, gate.duration]; kwargs...)
    solution.states[end]
end

function _objective_metrics(model, propagator, objective::UnitaryObjective, duration)
    basis = _logical_basis(model, objective.states)
    logical = adjoint(basis) * Matrix(propagator.data) * basis
    correction = _frame_correction(model, objective, duration)
    correction === nothing || (logical = correction * logical)
    target = objective.target
    dimension = size(target, 1)
    fidelity = abs(tr(adjoint(target) * logical))^2 / dimension^2
    survival = real(tr(adjoint(logical) * logical)) / dimension
    leakage = max(0.0, 1.0 - survival)
    loss = 1.0 - fidelity + objective.leakage_weight * leakage
    (;
        loss = Float64(real(loss)),
        fidelity = Float64(real(fidelity)),
        leakage = Float64(real(leakage)),
    )
end

function _objective_metrics(model, propagator, objective::StateTransferObjective, duration)
    initial = _objective_state(model, objective.initial)
    target = _objective_state(model, objective.target)
    final = Matrix(propagator.data) * initial
    fidelity = abs(dot(target, final))^2
    leakage = if objective.states === nothing
        0.0
    else
        basis = _logical_basis(model, objective.states)
        max(0.0, 1.0 - real(sum(abs2, adjoint(basis) * final)))
    end
    (;
        loss = Float64(real(1.0 - fidelity + leakage)),
        fidelity = Float64(real(fidelity)),
        leakage = Float64(real(leakage)),
    )
end

function _frame_correction(model, objective::UnitaryObjective, duration)
    frequencies = objective.frame_frequencies
    frequencies === nothing && return nothing
    subsystem_names = propertynames(model.spec.subsystems)
    propertynames(frequencies) == subsystem_names || error(
        "frame_frequencies keys $(propertynames(frequencies)) must match model " *
        "subsystems $subsystem_names in order.",
    )
    phases = ComplexF64[]
    frequency_values = Tuple(values(frequencies))
    for label in objective.states
        label isa Tuple && length(label) == length(frequency_values) || error(
            "Frame-corrected objective states must be full binary subsystem labels.",
        )
        all(level -> level in (0, 1), label) || error(
            "Frame-corrected objective states must contain only qubit levels 0 and 1.",
        )
        phase = π * duration * sum(
            frequency * (2level - 1) for (frequency, level) in zip(frequency_values, label)
        )
        push!(phases, cis(phase))
    end
    Diagonal(phases)
end

function _lab_frame_target(model, objective::UnitaryObjective, duration)
    correction = _frame_correction(model, objective, duration)
    correction === nothing ? objective.target : adjoint(correction) * objective.target
end

_objective_matrix(value::QuantumObject) = ComplexF64.(Matrix(value.data))
_objective_matrix(value::AbstractMatrix) = ComplexF64.(Matrix(value))
_objective_matrix(value) = error("Gate objective target must be a matrix or QuantumObject.")

function _logical_basis(model, labels)
    columns = [_objective_state(model, label) for label in labels]
    length(unique(labels)) == length(labels) ||
        error("Gate objective state labels must be unique.")
    basis = hcat(columns...)
    gram = adjoint(basis) * basis
    isapprox(gram, I; atol = 1e-8, rtol = 1e-8) ||
        error("Gate objective states must form an orthonormal basis.")
    basis
end

function _objective_state(model, state)
    vector = if state isa Tuple
        _labeled_model_state(model, state)
    elseif state isa QuantumObject
        vec(ComplexF64.(Array(state.data)))
    elseif state isa AbstractVector
        ComplexF64.(collect(state))
    else
        error("Gate objective states must be model labels, vectors, or QuantumObjects.")
    end
    length(vector) == size(model.hamiltonian, 1) || error(
        "Gate objective state has length $(length(vector)); expected " *
        "$(size(model.hamiltonian, 1)).",
    )
    norm(vector) > 0 || error("Gate objective states must be nonzero.")
    vector / norm(vector)
end

function _labeled_model_state(model, label)
    haskey(model.states, label) ||
        error("Model $(model.spec.name) does not contain state label $label.")
    if model.spec.states_to_keep === nothing
        return vec(ComplexF64.(Array(model.states[label].data)))
    end
    index = findfirst(==(label), model.spec.states_to_keep)
    index === nothing &&
        error("Model $(model.spec.name) does not retain state label $label.")
    vector = zeros(ComplexF64, length(model.spec.states_to_keep))
    vector[index] = 1
    vector
end
