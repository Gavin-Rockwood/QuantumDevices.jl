import OptimizationNLopt
import SciMLBase

"""Physical pulse fields keyed by `control_name/pulse/path`; derived durations are excluded."""
function calibration_parameters(gate::GateRecipe)
    result = Dict{String,Tuple}()
    for c in gate.controls, (path, (_, value)) in calibration_parameters(c.pulse)
        name = "$(c.name)/$path"
        _derived_pulse_path(c.pulse,split(path,'/')) && continue
        result[name] = (name,value)
    end
    return result
end

_matrix(x::AbstractMatrix) = Matrix(x)
_matrix(x::QuantumObject) = Matrix(x.data)
function _calibration_target(target)
    matrix = _matrix(target)
    n,m = size(matrix)
    n == m && n > 0 || throw(DimensionMismatch("target must be nonempty and square"))
    all(isfinite,matrix) && isapprox(matrix'*matrix,Matrix{ComplexF64}(I,n,n);atol=1e-8,rtol=1e-8) ||
        throw(ArgumentError("target must be unitary"))
    return matrix
end

"""
    gate_infidelity(target, actual)

One minus average pure-state target overlap, including loss from the selected
subspace: `1 - (tr(M'M) + abs2(tr(T'M)))/(d*(d+1))`. The target T is unitary;
M may leak and is never renormalized. Relative phases count, global phase does not.
"""
function gate_infidelity(target, actual)
    T = _calibration_target(target)
    M = _matrix(actual)
    size(T) == size(M) || throw(DimensionMismatch("target and gate matrix differ"))
    d = size(T,1)
    return clamp(real(1-(sum(abs2,M)+abs2(dot(T,M)))/(d*(d+1))),0.,1.)
end

struct CalibrationProblem{P,R,S,T,O,E,B,D,K}
    problem::P
    prepared::R
    coordinates::Tuple
    paths::Tuple{Vararg{String}}
    parameterization::NamedTuple
    initial_values::Tuple
    states::S
    target::T
    objective::O
    evolution_settings::E
    bounds::B
    derivative_settings::D
    state_spec::K
end

include("parameterization.jl")

# Resolve every output before reconstructing any physical pulse node.
function _candidate_recipe(start::GateRecipe,coordinates::Tuple,parameterization::NamedTuple,values)
    return _update_gate_pulses(start,_parameter_updates(coordinates,parameterization,values))
end

# Preserve selection order and distinguish default labels from explicit ket inputs.
function _capture_calibration_states(reference, requested)
    selections = requested === nothing ? deepcopy(_reference_states(reference)) : deepcopy(requested)
    _gate_kets(reference,selections)
    keys = if requested !== nothing
        all(state -> state isa Tuple,selections) ? Tuple(selections) : nothing
    elseif reference.tracking.tracked
        labels = Tuple(reference.tracking.keys)
        length(labels) == length(selections) &&
            all(index -> tracking_row(reference.tracking,labels[index]) == index,eachindex(labels)) ? labels : nothing
    else
        nothing
    end
    return selections,(;state_keys=keys,original_states=requested === nothing ? nothing : selections)
end

# Keep callback identities while independently capturing native backend/prototype data.
function capture_derivative_settings(adtype, function_kwargs)
    captured = map(value -> value isa Function ? value : deepcopy(value), (;function_kwargs...))
    return (;adtype=deepcopy(adtype), function_kwargs=captured)
end

"""
    CalibrationProblem(prepared, names, target; states=nothing,
                       objective=gate_infidelity, adtype=SciMLBase.NoAD(), function_kwargs=(;),
                       evolution_settings=(;), lb=nothing, ub=nothing, kwargs...)

Capture a starting gate and fixed reference projections. Direct names come from
`calibration_parameters(prepared.recipe)`. Alternatively use `CalibrationProblem(prepared,target;
variables=[Num=>initial,...],links=Dict(path=>expression,...))`: ordered coordinates
and bounds are separate from the affected physical paths. Real numeric variables are supported;
complex pulse values must be parameterized through real controls. `states` is an
ordered label/ket vector, defaulting to the full retained basis. The default
objective evolves all retained rows and projects only at measurement time, retaining
all physical phases. No reference-Hamiltonian inverse is applied.
A custom objective has contract `objective(target_matrix, projected_matrix)`.
`evolution_settings` forwards to `evolve`, including `remove_vacuum=true` for a
constant reference-vacuum energy origin (default false). This setting is captured
for both candidate and accepted-result reevaluation, not passed to either solver.
`adtype` and native `function_kwargs` (e.g. `grad`, `hess`, `cons`) construct the
SciML OptimizationFunction. Remaining keywords construct OptimizationProblem.
`solve` forwards the chosen backend and solve options unchanged, returning a CalibrationResult containing its
untouched native solution. `solve(problem.problem, ...)` is the raw escape. AD configuration does not guarantee end-to-end AD support.
The raw objective permits physically valid probes outside the optimization box;
`calibrated_recipe` and acceptance still enforce bounds. No clipping is applied.
Changing the reference recipe requires a fresh build/problem.
"""
# Rebase pulse definitions onto an existing endpoint without rebuilding its basis.
function _calibration_model(prepared::PreparedGate, model::NumericalModel)
    definition = draft(model.recipe.definition)
    merge!(definition.parameters,model.parameters)
    old = prepared.recipe.model_recipe
    old_scope = old isa CompositeRecipe ? old.scope : nothing
    new_scope = definition isa CompositeRecipe ? definition.scope : nothing
    terms = [copy_interaction(term,false;old_scope,new_scope) for term in prepared.recipe.terms]
    gate = capture_gate(nothing,GateRecipe(definition;controls=prepared.recipe.controls,
        terms,duration=prepared.recipe.duration_policy === :auto ? nothing : prepared.recipe.duration,
        metadata=prepared.recipe.metadata))
    extra,values = _reuse_gate_data(model,gate)
    expressions = vcat([model.hamiltonian],extra)
    required = union((refs(expression,:param) for expression in expressions)...)
    all(control -> control.target in required,gate.controls) ||
        throw(ArgumentError("control target must appear in replacement model Hamiltonian or gate terms"))
    foreach(control -> _resolve_pulse(control.pulse,values),gate.controls)
    for expression in expressions, (word,_) in physical_terms(expression)
        haskey(model.operator_products,word) || throw(ArgumentError(
            "complete operator product $word is not cached in replacement model; rebuild that model with the gate terms"))
    end
    return PreparedGate(gate,model,extra,values,nothing)
end

function CalibrationProblem(prepared::PreparedGate,paths,target;variables=missing,links=missing,kwargs...)
    variables === missing && links === missing || throw(ArgumentError("direct pulse paths cannot be combined with variables/links; use CalibrationProblem(prepared,target;variables,links)"))
    return _calibration_problem(prepared,paths,target;symbolic=false,kwargs...)
end

function CalibrationProblem(prepared::PreparedGate,target;variables=missing,links=missing,kwargs...)
    variables !== missing && links !== missing || throw(ArgumentError("symbolic calibration requires both variables and links"))
    return _calibration_problem(prepared,variables,target;symbolic=true,links,kwargs...)
end

function _calibration_problem(prepared::PreparedGate,selection,target;symbolic::Bool,links=nothing,
    states=missing,model::Union{Nothing,NumericalModel}=nothing,
    objective=gate_infidelity,adtype=SciMLBase.NoAD(),function_kwargs=(;),
    evolution_settings=(;),lb=nothing,ub=nothing,kwargs...)
    haskey(kwargs,:remove_idle) && throw(ArgumentError("unsupported keyword argument: remove_idle"))
    haskey(kwargs,:remove_vacuum) && throw(ArgumentError("pass remove_vacuum in evolution_settings, not as an optimizer keyword"))
    if model !== nothing
        states === missing && throw(ArgumentError("model replacement requires explicit states; coordinate kets cannot be transferred automatically"))
        prepared = _calibration_model(prepared,model)
    end
    states = states === missing ? nothing : states
    captured = deepcopy(prepared)
    normalized = symbolic ? _symbolic_parameterization(captured,selection,links) : _direct_parameterization(captured.recipe,selection)
    coordinates,paths,parameterization = normalized.coordinates,normalized.paths,normalized.parameterization
    initial = Float64.(collect(normalized.initial_values))
    lower = lb === nothing ? fill(-Inf,length(coordinates)) : Float64.(collect(lb))
    upper = ub === nothing ? fill(Inf,length(coordinates)) : Float64.(collect(ub))
    length(lower)==length(upper)==length(coordinates) || throw(DimensionMismatch("bounds must match selected parameters"))
    all(i -> !isnan(lower[i]) && !isnan(upper[i]) && lower[i]<=initial[i]<=upper[i],eachindex(initial)) ||
        throw(ArgumentError("invalid bounds or initial value outside bounds"))
    initial_recipe = _candidate_recipe(captured.recipe,coordinates,parameterization,initial)
    captured = PreparedGate(initial_recipe,captured.model,captured.terms,captured.parameters,captured.context)
    reference = _gate_reference(captured.model)
    selections,state_spec = _capture_calibration_states(reference,states)
    T = _calibration_target(target)
    size(T,1)==length(selections) || throw(DimensionMismatch("target does not match ordered states"))
    settings = deepcopy(merge((;remove_vacuum=false), (;evolution_settings...)))
    haskey(settings,:states) && throw(ArgumentError("pass states at CalibrationProblem construction, not inside evolution_settings"))
    settings.remove_vacuum isa Bool || throw(ArgumentError("evolution_settings.remove_vacuum must be a Bool"))
    derivative_settings = capture_derivative_settings(adtype,function_kwargs)
    rebuild = values -> _candidate_recipe(captured.recipe,coordinates,parameterization,values)
    loss = function(values,_)
        length(values)==length(coordinates) || throw(DimensionMismatch("wrong parameter count"))
        recipe = rebuild(values)
        candidate = PreparedGate(recipe,captured.model,captured.terms,captured.parameters,captured.context)
        result = evolve(candidate;states=selections,settings...)
        value = objective(T,gate_matrix(result))
        value isa Real && isfinite(value) || throw(ArgumentError("objective must return a finite real scalar"))
        return value
    end
    native_function = SciMLBase.OptimizationFunction(loss,derivative_settings.adtype;
        derivative_settings.function_kwargs...)
    bounded = lb !== nothing || ub !== nothing
    native = SciMLBase.OptimizationProblem(native_function,initial;
        lb=bounded ? lower : nothing,ub=bounded ? upper : nothing,kwargs...)
    return CalibrationProblem(native,captured,coordinates,paths,parameterization,Tuple(initial),selections,T,objective,
        settings,(;lb=Tuple(lower),ub=Tuple(upper)),derivative_settings,state_spec)
end

"""
    calibrated_recipe(problem::CalibrationProblem, values)

Materialize a concrete GateRecipe from coordinates in optimizer order. Resolve all
links and rebuild affected pulse fields atomically; do not optimize, evolve or mutate
problem snapshots. Invalid coordinates or pulse geometry throw before returning.
"""
function calibrated_recipe(p::CalibrationProblem,values)
    length(values)==length(p.coordinates) || throw(DimensionMismatch("wrong parameter count"))
    all(i -> p.bounds.lb[i]<=values[i]<=p.bounds.ub[i],eachindex(values)) || throw(ArgumentError("candidate outside bounds"))
    return _candidate_recipe(p.prepared.recipe,p.coordinates,p.parameterization,values)
end


include("CalibrationResult.jl")
include("QuantumGate.jl")
include("display.jl")
