"""
A candidate recipe and ordered target basis captured from a calibration solve.
`solution` is the untouched native optimizer result. Ordered `coordinates`, owned
`fitted_values` and `parameterization` retain the optimization family; `paths` lists
actual written physical fields and need not have the coordinate count. `validation` is nothing until explicitly reevaluated. Arrays remain
ordinary owned data; user callables may retain external mutable state.
"""
struct CalibrationResult{R,S,K,U,T,V,C,P,E}
    recipe::R
    solution::S
    coordinates::Tuple
    paths::Tuple{Vararg{String}}
    parameterization::NamedTuple
    fitted_values::Tuple
    state_keys::K
    states::U
    target::T
    virtual_z::V
    settings::C
    provenance::P
    validation::E
end

# Capture native problem fields once; solver-specific result data stays in solution.
function _calibration_result_settings(problem::CalibrationProblem)
    native = problem.problem
    options = map(value -> value isa Function ? value : deepcopy(value),(;native.kwargs...))
    return (;evolution=deepcopy(problem.evolution_settings),
        derivatives=capture_derivative_settings(problem.derivative_settings.adtype,
                                               problem.derivative_settings.function_kwargs),
        problem=(;initial_values=problem.initial_values,lb=deepcopy(native.lb),ub=deepcopy(native.ub),
            int=deepcopy(native.int),lcons=deepcopy(native.lcons),ucons=deepcopy(native.ucons),
            sense=native.sense,kwargs=options))
end

"""
    CalibrationResult(problem, solution)

Capture a feasible pulse candidate without evolving or rebuilding its model.
Keep `solution` by identity, including unsuccessful return codes. Unusable or
out-of-box parameters throw; obtain raw optimizer results with
`solve(problem.problem, algorithm; ...)`. No virtual-Z correction is fitted.
"""
function CalibrationResult(problem::CalibrationProblem, solution)
    hasproperty(solution,:u) && hasproperty(solution,:retcode) ||
        throw(ArgumentError("solution must expose u and retcode"))
    fitted_values = Tuple(deepcopy(collect(solution.u)))
    recipe = capture_gate(nothing,calibrated_recipe(problem,fitted_values))
    source = problem.prepared
    prepared = PreparedGate(recipe,deepcopy(source.model),deepcopy(source.terms),
                            copy(source.parameters),deepcopy(source.context))
    reference = _gate_reference(prepared.model;parameters=prepared.parameters)
    states = _gate_kets(reference,problem.states)
    provenance = (;prepared,objective=problem.objective,
        original_states=deepcopy(problem.state_spec.original_states))
    return CalibrationResult(recipe,solution,deepcopy(problem.coordinates),problem.paths,
        deepcopy(problem.parameterization),fitted_values,
        deepcopy(problem.state_spec.state_keys),states,copy(problem.target),nothing,
        _calibration_result_settings(problem),provenance,nothing)
end

"""
    solve(problem::CalibrationProblem[, algorithm]; kwargs...)

Return a CalibrationResult with the untouched native optimizer solution. Without
an algorithm, use NLopt's bound-aware LN_NELDERMEAD; its `maxiters` limits objective
evaluations. Explicit algorithms retain their native SciML keyword semantics.
No additional evolution or coordinate transformation is performed.
"""
function SciMLBase.solve(problem::CalibrationProblem; kwargs...)
    return SciMLBase.solve(problem,OptimizationNLopt.NLopt.LN_NELDERMEAD;kwargs...)
end

# Only this wrapper changes return packaging; the underlying native problem stays raw.
function SciMLBase.solve(problem::CalibrationProblem, args...; kwargs...)
    solution = SciMLBase.solve(problem.problem,args...;kwargs...)
    return CalibrationResult(problem,solution)
end

"""
    CalibrationProblem(result; kwargs...)
    CalibrationProblem(result, paths, target=result.target; kwargs...)

Warm-start from the captured candidate without rebuilding its physical model.
Evolution settings merge with the original solve settings; explicit
`function_kwargs` replace derivative hooks. Explicit paths or mapping declarations require explicit
`lb`, `ub`, and `function_kwargs`, and reset integrality/constraint bounds.
"""
function CalibrationProblem(result::CalibrationResult;target=result.target,variables=missing,links=missing,kwargs...)
    reset = variables !== missing || links !== missing
    if reset
        variables !== missing && links !== missing || throw(ArgumentError("replace symbolic coordinates with both variables and links"))
        return _continue_calibration(result,variables,target;symbolic=true,links,reset=true,kwargs...)
    end
    symbolic = result.parameterization.mode === :symbolic
    selection = result.coordinates
    # Construct pairs explicitly; never infer an inverse from physical pulse values.
    symbolic && (selection = [variable=>value for (variable,value) in zip(result.coordinates,result.fitted_values)])
    return _continue_calibration(result,selection,target;symbolic,
        links=symbolic ? result.parameterization.links : nothing,reset=false,kwargs...)
end

function CalibrationProblem(result::CalibrationResult,paths,target=result.target;variables=missing,links=missing,kwargs...)
    variables === missing && links === missing || throw(ArgumentError("direct path selection cannot be combined with variables/links"))
    return _continue_calibration(result,paths,target;symbolic=false,reset=true,kwargs...)
end

function _continue_calibration(result,selection,target;symbolic::Bool,links=nothing,reset::Bool,kwargs...)
    options = (;kwargs...)
    if get(options,:model,nothing) !== nothing && isnothing(result.state_keys) && !haskey(options,:states)
        throw(ArgumentError("replacement model requires explicit states for a result without occupation labels"))
    end
    if reset
        all(key -> haskey(options,key),(:lb,:ub,:function_kwargs)) || throw(ArgumentError(
            "explicit coordinates/mapping require lb, ub, and function_kwargs; use nothing, nothing, and (;) to reset"))
    end
    derivatives = result.settings.derivatives
    if haskey(options,:adtype) && !isempty(derivatives.function_kwargs) && !haskey(options,:function_kwargs)
        throw(ArgumentError("an explicit adtype with inherited derivative hooks requires explicit function_kwargs"))
    end
    native = result.settings.problem
    indexed = reset ? (;int=nothing,lcons=nothing,ucons=nothing) :
        (;int=deepcopy(native.int),lcons=deepcopy(native.lcons),ucons=deepcopy(native.ucons))
    selected_states = isnothing(result.state_keys) ? deepcopy(result.states) : collect(result.state_keys)
    defaults = merge(native.kwargs,indexed,(;sense=native.sense,
        lb=deepcopy(native.lb),ub=deepcopy(native.ub),states=selected_states,
        objective=result.provenance.objective,adtype=derivatives.adtype,
        function_kwargs=derivatives.function_kwargs))
    evolution = merge(result.settings.evolution,get(options,:evolution_settings,(;)))
    arguments = merge(defaults,options,(;evolution_settings=evolution))
    source = result.provenance.prepared
    # Restore written fields from owned values even when native u has been mutated.
    recipe = _candidate_recipe(result.recipe,result.coordinates,result.parameterization,result.fitted_values)
    prepared = PreparedGate(recipe,source.model,source.terms,source.parameters,source.context)
    return _calibration_problem(prepared,selection,target;symbolic,links,arguments...)
end

# Replay the fixed snapshot; explicit caller options override captured evolution settings.
function evolve(result::CalibrationResult; kwargs...)
    options = merge((;states=result.states),result.settings.evolution,(;kwargs...))
    return evolve(result.provenance.prepared;options...)
end

# A validated result owns its data independently, while preserving native solution identity.
function _validated_calibration_result(result, validation)
    data = deepcopy((;recipe=result.recipe,coordinates=result.coordinates,paths=result.paths,
        parameterization=result.parameterization,fitted_values=result.fitted_values,
        state_keys=result.state_keys,states=result.states,
        target=result.target,virtual_z=result.virtual_z,prepared=result.provenance.prepared,
        original_states=result.provenance.original_states,validation))
    derivatives = capture_derivative_settings(result.settings.derivatives.adtype,
                                              result.settings.derivatives.function_kwargs)
    problem = deepcopy(result.settings.problem)
    options = map(value -> value isa Function ? value : deepcopy(value),result.settings.problem.kwargs)
    settings = (;evolution=deepcopy(result.settings.evolution),derivatives,
                 problem=merge(problem,(;kwargs=options)))
    provenance = (;prepared=data.prepared,objective=result.provenance.objective,
                   original_states=data.original_states)
    return CalibrationResult(data.recipe,result.solution,data.coordinates,data.paths,
        data.parameterization,data.fitted_values,data.state_keys,
        data.states,data.target,data.virtual_z,settings,provenance,data.validation)
end

"""
    validate_calibration(result; evolution_settings=(;), save_evolution=false)

Perform one explicit reevaluation with the captured target/input order. Return a
new result with `validation.loss`, its effective settings, and optional saved
evolution. The optimizer-reported objective and original result are unchanged.
"""
function validate_calibration(result::CalibrationResult; evolution_settings=(;), save_evolution=false)
    haskey(evolution_settings,:states) && throw(ArgumentError("validation keeps the captured target state order"))
    settings = merge(result.settings.evolution,(;evolution_settings...))
    evolution = evolve(result;settings...)
    loss = result.provenance.objective(result.target,gate_matrix(evolution))
    loss isa Real && isfinite(loss) || throw(ArgumentError("validation objective must be finite and real"))
    validation = (;loss,evolution_settings=deepcopy(settings),evolution=save_evolution ? evolution : nothing)
    return _validated_calibration_result(result,validation)
end

