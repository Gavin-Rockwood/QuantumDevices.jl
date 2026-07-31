module QuantumDevicesPiccoloExt

using LinearAlgebra
using SparseArrays
import Piccolo
import QuantumDevices

const QD = QuantumDevices

function QD.optimize_gate(
    problem::QD.PiccoloOptimizationProblem;
    template_kwargs = Dict{Symbol,Any}(),
    solve_kwargs = Dict{Symbol,Any}(),
)
    _run_piccolo(problem; template_kwargs, solve_kwargs)
end

function QD._optimize_piccolo(
    gate::QD.GateSpec,
    objective::QD.AbstractGateObjective,
    selectors,
    bounds;
    solver_kwargs = Dict{Symbol,Any}(),
    knots::Integer = 100,
    final_fidelity::Real = 0.999,
    template_kwargs = Dict{Symbol,Any}(),
    solve_kwargs = Dict{Symbol,Any}(),
)
    path_selectors = QD.ParamPath[]
    optimize_duration = false
    for selector in selectors
        if selector === :duration
            optimize_duration = true
        elseif selector isa QD.ParamPath
            push!(path_selectors, selector)
        elseif selector isa Symbol
            error(
                "Piccolo accepts model ParamPath controls and optional :duration; " *
                "named recipe parameter $selector is unsupported.",
            )
        end
    end
    isempty(path_selectors) &&
        error("Piccolo optimization requires at least one ParamPath control.")
    0 < final_fidelity <= 1 || error("final_fidelity must lie in (0, 1].")
    if optimize_duration
        objective isa QD.UnitaryObjective &&
            objective.frame_frequencies !== nothing &&
            error(
                "Piccolo minimum-time optimization is unavailable for frame-corrected " *
                "unitary targets because the lab-frame target depends on duration.",
            )
        lower, upper = bounds[:duration]
        lower <= gate.duration <= upper ||
            error("The gate duration must lie within the Piccolo duration bounds.")
    end
    selected = Set(path_selectors)
    for (path, control) in gate.controls
        path in selected && continue
        default = QD._model_value(gate.modelspec, path)
        nondefault =
            applicable(control, 0.0) ||
            !isapprox(control, QD._at(default, 0.0); atol = 1e-10, rtol = 1e-10)
        nondefault && error(
            "Piccolo control $(join(path.parts, "/")) is time-dependent and non-default " *
            "but was not selected for optimization.",
        )
    end
    path_bounds = Dict(path => bounds[path] for path in path_selectors)
    ordered_paths = sort(copy(path_selectors); by = path -> path.parts)
    initial = Matrix{Float64}(undef, length(ordered_paths), Int(knots))
    control_modulations = Dict{QD.ParamPath,Any}()
    times = range(0.0, gate.duration; length = Int(knots))
    for (index, path) in enumerate(ordered_paths)
        carrier = QD._carrier_envelope(gate, path)
        if carrier === nothing
            control = get(gate.controls, path, QD._model_value(gate.modelspec, path))
            initial[index, :] .= Float64[QD._at(control, time) for time in times]
        else
            initial[index, :] .= Float64[carrier.envelope(time) for time in times]
            control_modulations[path] = carrier.modulation
        end
    end
    problem = QD.PiccoloOptimizationProblem(
        gate.name,
        gate.modelspec,
        path_bounds,
        objective;
        duration = gate.duration,
        knots,
        initial_controls = initial,
        control_modulations,
        template_kwargs,
        solve_kwargs,
    )
    duration_bounds = optimize_duration ? bounds[:duration] : nothing
    result = _run_piccolo(
        problem;
        duration_bounds,
        final_fidelity = Float64(final_fidelity),
        solver_kwargs,
    )
    if optimize_duration
        result.minimizer[:duration] = result.gate.duration
    end
    result
end

function _run_piccolo(
    problem::QD.PiccoloOptimizationProblem;
    template_kwargs = Dict{Symbol,Any}(),
    solve_kwargs = Dict{Symbol,Any}(),
    duration_bounds = nothing,
    final_fidelity::Float64 = 0.999,
    solver_kwargs = Dict{Symbol,Any}(),
)
    built = QD.model(problem.modelspec)
    decomposition = QD._model_control_decomposition(built, problem.controls)
    drives = Piccolo.AbstractDrive[]
    for term in decomposition.terms
        operator = sparse(2π .* term.operator)
        if term.linear_weights === nothing
            modulated_paths = problem.controls[
                findall(modulation -> modulation !== nothing, problem.control_modulations)
            ]
            isempty(modulated_paths) || error(
                "Piccolo cannot apply carrier modulation to nonlinear control term(s) " *
                "involving $(join(string.(modulated_paths), ", ")).",
            )
            push!(
                drives,
                Piccolo.NonlinearDrive(
                    operator,
                    term.coefficient;
                    active_controls = collect(eachindex(problem.controls)),
                ),
            )
        else
            for (index, weight) in enumerate(term.linear_weights)
                iszero(weight) && continue
                drive = Piccolo.LinearDrive(weight .* operator, index)
                modulation = problem.control_modulations[index]
                push!(
                    drives,
                    modulation === nothing ? drive : Piccolo.ModulatedDrive(drive, modulation),
                )
            end
        end
    end
    isempty(drives) &&
        error("Piccolo conversion found no control-dependent Hamiltonian terms.")
    system =
        Piccolo.QuantumSystem(sparse(2π .* decomposition.drift), drives, problem.bounds)

    times = collect(range(0.0, problem.duration; length = problem.knots))
    pulse = Piccolo.ZeroOrderPulse(
        copy(problem.initial_controls),
        times;
        initial_value = copy(problem.initial_controls[:, 1]),
        final_value = copy(problem.initial_controls[:, end]),
    )
    trajectory = _piccolo_trajectory(system, pulse, built, problem.objective)
    step = problem.duration / (problem.knots - 1)
    options = Piccolo.PiccoloOptions(timesteps_all_equal = true, display = :silent)
    template = Dict{Symbol,Any}(
        :Q => 100.0,
        :R => 1e-2,
        :Δt_bounds => (step, step),
        :piccolo_options => options,
    )
    merge!(template, problem.template_kwargs)
    merge!(template, template_kwargs)
    template[:Δt_bounds] = (step, step)
    qcp = Piccolo.SmoothPulseProblem(trajectory, problem.knots; template...)

    solve_options = copy(problem.solve_kwargs)
    merge!(solve_options, solve_kwargs)
    get!(solve_options, :max_iter, 100)
    get!(solve_options, :verbose, false)
    get!(solve_options, :print_level, 0)
    Piccolo.solve!(qcp; solve_options...)

    template_name = :SmoothPulseProblem
    if duration_bounds !== nothing
        problem.objective isa QD.UnitaryObjective &&
            problem.objective.frame_frequencies !== nothing &&
            error(
                "Piccolo minimum-time optimization is unavailable for frame-corrected " *
                "unitary targets.",
            )
        lower, upper = duration_bounds
        step_bounds = (lower / (problem.knots - 1), upper / (problem.knots - 1))
        qcp = Piccolo.MinimumTimeProblem(
            qcp;
            final_fidelity,
            Δt_bounds = step_bounds,
            piccolo_options = options,
        )
        Piccolo.solve!(qcp; solve_options...)
        template_name = :MinimumTimeProblem
    end

    optimized_pulse = Piccolo.get_pulse(qcp.qtraj)
    optimized_duration = Float64(Piccolo.duration(qcp.qtraj))
    pulse_samples, _ = Piccolo.sample(optimized_pulse, problem.knots)
    controls = Dict{QD.ParamPath,Any}()
    for (index, path) in enumerate(problem.controls)
        default = problem.initial_controls[index, 1]
        modulation = problem.control_modulations[index]
        pulse_samples[index, 1] = default
        pulse_samples[index, end] = default
        controls[path] =
            let samples = copy(pulse_samples[index, :]),
                duration = optimized_duration,
                knots = problem.knots,
                default = default

                time -> begin
                    time <= 0 && return default
                    time >= duration && return default
                    knot =
                        clamp(floor(Int, time * (knots - 1) / duration) + 1, 1, knots)
                    value = samples[knot]
                    modulation === nothing ? value : value * modulation(time)
                end
            end
    end
    gate = QD.GateSpec(
        problem.name,
        problem.modelspec;
        duration = optimized_duration,
        controls,
    )
    metrics = QD._evaluate_gate_objective(gate, problem.objective; solver_kwargs)
    settings = Dict{Symbol,Any}(
        :backend => :Piccolo,
        :template => template_name,
        :knots => problem.knots,
        :control_paths => copy(problem.controls),
        :initial_controls => copy(problem.initial_controls),
        :carrier_paths => [
            path for (path, modulation) in
            zip(problem.controls, problem.control_modulations) if modulation !== nothing
        ],
        :template_kwargs => template,
        :solve_kwargs => solve_options,
        :duration_bounds => duration_bounds,
        :final_fidelity => final_fidelity,
    )
    QD.GateOptimizationResult(
        gate,
        Dict{Symbol,Float64}(),
        metrics.loss,
        metrics.fidelity,
        metrics.leakage,
        true,
        0,
        settings,
        NamedTuple[],
        qcp,
    )
end

function _piccolo_trajectory(system, pulse, model, objective::QD.UnitaryObjective)
    basis = QD._logical_basis(model, objective.states)
    dimension = size(model.hamiltonian, 1)
    projector = basis * adjoint(basis)
    logical_target = QD._lab_frame_target(model, objective, Piccolo.duration(pulse))
    target =
        basis * logical_target * adjoint(basis) +
        Matrix{ComplexF64}(I, dimension, dimension) - projector
    isapprox(target * adjoint(target), I; atol = 1e-8, rtol = 1e-8) ||
        error("Piccolo unitary target is not unitary on the selected model space.")
    Piccolo.UnitaryTrajectory(system, pulse, target)
end

function _piccolo_trajectory(system, pulse, model, objective::QD.StateTransferObjective)
    initial = QD._objective_state(model, objective.initial)
    target = QD._objective_state(model, objective.target)
    Piccolo.KetTrajectory(system, pulse, initial, target)
end

end
