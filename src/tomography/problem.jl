"""
    ModelTomographyProblem(model, paths; predict, observed, uncertainty=1, loss=:auto, kwargs...)
    ModelTomographyProblem(model, paths; sequences, preparation, measurement, observed=nothing, counts=nothing, shots=nothing, kwargs...)
    ModelTomographyProblem(model, data, paths; loss=:auto,
        probability_floor=1e-12, evolution_kwargs=(; abstol=1e-8, reltol=1e-8), threaded::Bool=false, subspace=nothing, leakage_outcome=nothing,
    normalize::Bool=true, parameter_scales=nothing, scale_function=nothing, lb=nothing, ub=nothing, threaded_components::Bool=true, normalize_probabilities::Bool=false, reference_model=:initial, tracking_steps=0, kwargs...)

For sequence fits, pass plain `sequences` and shared `preparation`/`measurement`
with exactly one of `observed` probability vectors or `counts` vectors, ordered by
sequence and measurement outcome. Optional `shots` is a shared integer or a
vector in sequence order (entries may be `nothing`). Counts infer their shot totals.
Do not mix scalar `predict` with `sequences`. Existing per-record inputs remain
available for experiments with different preparations or measurements.

Fit finite real model/gate scalars selected by writable slash paths. Accepts
[`SequenceData`](@ref) records, or a scalar `predict(model)` function and an
`observed` NamedTuple. Scalar predictions must return a NamedTuple with exactly
the observation names, in any order. Values must be finite real scalars.
`uncertainty` is a nonnegative scalar (default 1) or a NamedTuple with those names.
Zero uses an unweighted residual for that observation, rather than division by
zero or an exact equality constraint.
Residuals follow observation order. Exposes a
standard SciML optimization `.problem`, immutable `.rebuild(values)`, and `.paths`.
Solve using `solve(setup, algorithm; kwargs...)`; no optimizer is selected here.
By default `normalize=true` uses `physical = initial + scale .* coordinates`,
with zero optimizer initial coordinates. Scales are absolute initial values,
or 1 for zero-valued parameters. `parameter_scales` overrides scales in path
order, or `scale_function(path, initial_value)` supplies them. Scales stay fixed.
Bounds are supplied in physical units and transformed automatically. With
`normalize=false`, optimizer coordinates are physical values. Use
[`physical_parameters`](@ref) to interpret `solution.u`; `.rebuild` and
`tomography_residuals` consume optimizer coordinates in both modes.
Candidate values are applied together: each changed component is reconstructed
once. Physical candidates use `make_model`; no incremental model cache is kept.
`threaded_components=true` reconstructs changed components concurrently on Julia threads. Set it to `false` for serial
reconstruction; custom component constructors must be thread-safe.
Gate-only candidates reuse the idle model. `reference_model=:initial` and
`tracking_steps=0` match each candidate once against the fixed initial eigenbasis,
transporting rebuilt component bases for physical overlaps. Positive tracking
steps add intermediate models. Automatic initial references fall back to bare
tracking when no compatible basis transport exists; explicit references must
be compatible. Set `reference_model=nothing, tracking_steps=18` for the original
bare-state tracking. Physical rebuilds do not inherit gates.
Sequence fits with nonempty gates require an explicit `gates` dictionary. Selected
gate paths update that dictionary independently of physical reconstruction. Every candidate starts from the input
model, independently of earlier optimizer evaluations. Paths must resolve to
distinct writable destinations; aliases and overlapping paths are rejected.
`loss=:auto` selects `:likelihood` when every record has shots and `:least_squares`
when none does. Mixed shot-qualified sequence datasets require explicit
`:least_squares`. Scalar datasets always use least squares, with residuals
`(predicted - observed) / uncertainty` for positive uncertainty, or the raw
difference for zero uncertainty. They reconstruct each candidate
once and evaluate all predictions without gate evolution.
Likelihood is summed multinomial negative log likelihood without constants;
probabilities with shots use `-N*sum(f .* log(p))`, without rounding counts.
Only logarithms use `probability_floor`; zero-count terms contribute zero.
Least squares sums squared outcome-probability residuals without observable-
variance weighting. Additional keywords, including bounds, pass to SciML.
`threaded=true` builds distinct gate propagators in parallel on Julia threads;
it and `evolution_kwargs` have no effect on scalar fits.
`subspace` and `leakage_outcome` enable projected sequence fitting as described
by [`predict_sequence`](@ref). `normalize_probabilities=true` normalizes predicted
outcome weights after complementary leakage is assigned; the default `false`
retains strict validation. This option has no effect on scalar fits.
For sequence fits, `diagnostics=report -> ...` receives per-prediction wall-clock
timings in seconds: rebuilding, gate updates, input setup, propagation (including
per-gate pairs), and sequence/measurement evaluation, plus physical parameters.
Constructor validation is included; first calls may include compilation. The
callback runs after prediction and is excluded from the reported total.
Shared controls and solver callbacks must be safe for concurrent calls.
Automatic differentiation is not guaranteed through model reconstruction.
"""
struct ModelTomographyProblem{P,F,R,S,T}
    problem::P
    rebuild::F
    residuals::R
    paths::S
    rescale::T
end

function _tomography_parameter_setup(model, paths; normalize::Bool=true,
    parameter_scales=nothing, scale_function=nothing, lb=nothing, ub=nothing, threaded_components::Bool=true, reference_model=:initial, tracking_steps=0)
    paths isa AbstractString && throw(ArgumentError("Pass a collection of parameter paths"))
    paths = collect(paths)
    !isempty(paths) && all(p -> p isa AbstractString, paths) && allunique(paths) ||
        throw(ArgumentError("Parameter paths must be a nonempty unique collection of strings"))
    _tracking_step_count(tracking_steps,nothing)
    reference_model === :initial || reference_model === nothing || reference_model isa DeviceModel ||
        throw(ArgumentError("reference_model must be :initial, nothing, or a DeviceModel"))
    update_tree = _prepare_model_updates(model, paths)
    initial = Float64[]
    for path in paths
        value = getpath(model, path)
        _tomography_real_scalar(value) || throw(ArgumentError("Parameter $path must be finite and real"))
        push!(initial, value)
    end
    parameter_scales === nothing || scale_function === nothing ||
        throw(ArgumentError("Supply parameter_scales or scale_function, not both"))
    normalize || (parameter_scales === nothing && scale_function === nothing) ||
        throw(ArgumentError("Parameter scales require normalize=true"))
    scales = if !normalize
        ones(length(initial))
    elseif parameter_scales !== nothing
        supplied = collect(parameter_scales)
        supplied isa AbstractVector || throw(ArgumentError("parameter_scales must be a vector in path order"))
        length(supplied) == length(initial) || throw(DimensionMismatch("One scale is needed per path"))
        supplied
    elseif scale_function !== nothing
        [scale_function(path, value) for (path, value) in zip(paths, initial)]
    else
        [iszero(value) ? 1.0 : abs(value) for value in initial]
    end
    for (path, scale) in zip(paths, scales)
        _tomography_real_scalar(scale) && scale > 0 ||
            throw(ArgumentError("Scale for $path must be finite and positive"))
    end
    scales = Float64.(scales)
    all(s -> isfinite(s) && s > 0, scales) || throw(ArgumentError("Scales must fit positive Float64 values"))
    offset = normalize ? initial : zeros(length(initial))
    rescale = function(values)
        length(values) == length(paths) || throw(DimensionMismatch("Candidate length does not match selected paths"))
        all(_tomography_real_scalar, values) || throw(ArgumentError("Candidate values must be finite real scalars"))
        physical = offset .+ scales .* values
        all(_tomography_real_scalar, physical) || throw(ArgumentError("Rescaled parameters must be finite"))
        return physical
    end
    rebuild = function(values)
        return _apply_model_updates(model,update_tree,rescale(values);
            threaded_components,reference_model,tracking_steps)
    end
    function transform_bound(bound)
        bound === nothing && return nothing
        bound = collect(bound)
        bound isa AbstractVector || throw(ArgumentError("Bounds must be vectors in path order"))
        length(bound) == length(paths) || throw(DimensionMismatch("Bounds must match selected paths"))
        all(v -> v isa Real && !(v isa Bool) && !isnan(v), bound) ||
            throw(ArgumentError("Bounds must be real and not NaN"))
        return (bound .- offset) ./ scales
    end
    lower, upper = transform_bound(lb), transform_bound(ub)
    coordinates = normalize ? zeros(length(initial)) : copy(initial)
    for i in eachindex(initial)
        lo = lower === nothing ? -Inf : lower[i]
        hi = upper === nothing ? Inf : upper[i]
        lo <= coordinates[i] <= hi || throw(ArgumentError("Initial parameter $(paths[i]) must lie within bounds"))
    end
    return paths, coordinates, rebuild, rescale, (; lb=lower, ub=upper)

end

function ModelTomographyProblem(model::DeviceModel, paths; predict=nothing, observed=nothing,
    sequences=nothing, preparation=nothing, measurement=nothing, counts=nothing, shots=nothing, kwargs...)
    if sequences === nothing
        preparation === nothing && measurement === nothing && counts === nothing && shots === nothing ||
            throw(ArgumentError("Sequence inputs require sequences"))
        predict === nothing && throw(ArgumentError("Scalar fitting requires predict; sequence fitting requires sequences"))
        return _scalar_tomography_problem(model,paths;predict,observed,kwargs...)
    end
    predict === nothing || throw(ArgumentError("Supply predict or sequences, not both"))
    preparation !== nothing && measurement !== nothing ||
        throw(ArgumentError("Sequence fitting requires preparation and measurement"))
    data = _tomography_sequence_data(sequences,preparation,measurement;observed,counts,shots)
    return ModelTomographyProblem(model,data,paths;kwargs...)
end

function _scalar_tomography_problem(model::DeviceModel, paths; predict, observed, uncertainty=1,
    loss=:auto, threaded::Bool=false, evolution_kwargs=(;), normalize::Bool=true,
    parameter_scales=nothing, scale_function=nothing, lb=nothing, ub=nothing, threaded_components::Bool=true, normalize_probabilities::Bool=false, reference_model=:initial, tracking_steps=0, kwargs...)
    paths, initial, rebuild, rescale, bounds = _tomography_parameter_setup(model, paths;
        normalize, parameter_scales, scale_function, lb, ub, threaded_components, reference_model, tracking_steps)
    observed isa NamedTuple && !isempty(observed) ||
        throw(ArgumentError("Observed values must be a nonempty NamedTuple"))
    names = keys(observed)
    targets = _tomography_named_values(observed, names, "Observation")
    scales = if uncertainty isa NamedTuple
        _tomography_named_values(uncertainty, names, "Uncertainty")
    else
        _tomography_real_scalar(uncertainty) && uncertainty >= 0 ||
            throw(ArgumentError("Uncertainty for $(names) must be a finite nonnegative real scalar"))
        map(_ -> uncertainty, names)
    end
    for (name, scale) in zip(names, scales)
        scale >= 0 || throw(ArgumentError("Uncertainty for $name must be nonnegative"))
    end
    scales = map(scale -> iszero(scale) ? one(scale) : scale, scales)
    loss in (:auto, :least_squares) || throw(ArgumentError("Scalar fits require least squares"))
    residuals = function(values)
        predictions = _tomography_named_values(predict(rebuild(values)), names, "Prediction")
        return collect((predictions .- targets) ./ scales)
    end
    residuals(initial)
    objective = (values, _) -> sum(abs2, residuals(values))
    problem = SciMLBase.OptimizationProblem(objective, initial; bounds..., kwargs...)
    return ModelTomographyProblem(problem, rebuild, residuals, paths, rescale)
end

function ModelTomographyProblem(model::DeviceModel, data, paths; loss=:auto,
    probability_floor=1e-12, evolution_kwargs=(; abstol=1e-8, reltol=1e-8), threaded::Bool=false, subspace=nothing, leakage_outcome=nothing,
    normalize::Bool=true, parameter_scales=nothing, scale_function=nothing, lb=nothing, ub=nothing, threaded_components::Bool=true, normalize_probabilities::Bool=false, gates=nothing, reference_model=:initial, tracking_steps=0, diagnostics=nothing, kwargs...)
    data = collect(data)
    !isempty(data) || throw(ArgumentError("A tomography problem needs observations"))
    all(d -> d isa SequenceData,data) || throw(ArgumentError("Expected SequenceData records"))
    gates === nothing && any(d -> !isempty(d.sequence),data) &&
        throw(ArgumentError("Sequence fitting requires explicit gates=..."))
    gates = gates === nothing ? Dict{Any,Any}() : copy(gates)
    model = _copy_model_fields(model,Dict(:gates=>gates))
    paths, initial, rebuild, rescale, bounds = _tomography_parameter_setup(model, paths;
        normalize, parameter_scales, scale_function, lb, ub, threaded_components, reference_model, tracking_steps)
    data = collect(data)
    !isempty(data) || throw(ArgumentError("A tomography problem needs observations"))
    all(d -> d isa SequenceData, data) ||
        throw(ArgumentError("Expected SequenceData records; for scalar fits use ModelTomographyProblem(model, paths; predict, observed)"))
    qualified = [d.shots !== nothing for d in data]
    if loss === :auto
        loss = all(qualified) ? :likelihood : !any(qualified) ? :least_squares :
            throw(ArgumentError("Mixed shot-qualified data require loss=:least_squares"))
    end
    loss in (:likelihood, :least_squares) || throw(ArgumentError("Unknown tomography loss: $loss"))
    loss === :likelihood && !all(qualified) && throw(ArgumentError("Likelihood requires shots for every record"))
    isfinite(probability_floor) && 0 < probability_floor < 1 || throw(ArgumentError("Probability floor must lie strictly between zero and one"))
    gate_tree = get(_prepare_model_updates(model,paths).children,:gates,nothing)
    candidate_gates = values -> gate_tree === nothing ? gates :
        _apply_model_updates(gates,gate_tree,rescale(values))
    predictions = function(values)
        diagnostics === nothing && return predict_sequences(rebuild(values),data; evolution_kwargs,threaded,
            subspace,leakage_outcome,normalize_probabilities,gates=candidate_gates(values))
        started = time_ns()
        candidate = rebuild(values)
        rebuild_seconds = (time_ns()-started)/1e9
        gate_started = time_ns()
        selected_gates = candidate_gates(values)
        gate_update_seconds = (time_ns()-gate_started)/1e9
        return predict_sequences(candidate,data; evolution_kwargs,threaded,subspace,
            leakage_outcome,normalize_probabilities,gates=selected_gates,
            diagnostics=report -> diagnostics(merge(report,(; rebuild_seconds, gate_update_seconds,
                parameters=copy(rescale(values)), total_seconds=(time_ns()-started)/1e9))))
    end
    function checked_predictions(values)
        predicted = predictions(values)
        for (p, datum) in zip(predicted, data)
            length(p) == length(datum.probabilities) || throw(DimensionMismatch("Observed and predicted outcome counts differ"))
        end
        return predicted
    end
    checked_predictions(initial)
    residuals = values -> vcat((p .- d.probabilities for (p, d) in zip(checked_predictions(values), data))...)
    objective = if loss === :least_squares
        (values, _) -> sum(abs2, residuals(values))
    else
        function(values, _)
            predicted = checked_predictions(values)
            return sum(zip(predicted, data)) do (p, datum)
                weights = datum.counts === nothing ? datum.shots .* datum.probabilities : datum.counts
                sum(w == 0 ? zero(p[k]) : -w * log(max(p[k], probability_floor)) for (k, w) in enumerate(weights))
            end
        end
    end
    problem = SciMLBase.OptimizationProblem(objective, initial; bounds..., kwargs...)
    return ModelTomographyProblem(problem, rebuild, loss === :least_squares ? residuals : nothing, paths, rescale)
end

"""
    reconstructed_model(problem::ModelTomographyProblem, values)

Rebuild a fitted model from optimizer coordinates, e.g. `solution.u`, without
modifying the input. Normalized coordinates are rescaled automatically.
Inspect the optimizer retcode and fit quality before accepting its values.
"""
reconstructed_model(problem::ModelTomographyProblem, values) = problem.rebuild(values)

"""
    tomography_residuals(problem::ModelTomographyProblem, values)

Return prediction-minus-observation residuals in dataset and outcome order.
Inputs are optimizer coordinates, rescaled before predicting observations.
For scalar fits, residuals follow observation order and are divided by positive
uncertainty; zero uncertainty leaves the residual unweighted.
Only available for least squares; their squared norm equals the scalar objective.
External residual-based solvers may call this function directly.
"""
function tomography_residuals(problem::ModelTomographyProblem, values)
    problem.residuals === nothing && throw(ArgumentError("Likelihood problems do not have least-squares residuals"))
    return problem.residuals(values)
end

"""
    physical_parameters(problem::ModelTomographyProblem, values)

Convert optimizer coordinates (e.g. `solution.u`) to physical values in selected
path order. With `normalize=false`, return a copy of the supplied physical values.
"""
physical_parameters(problem::ModelTomographyProblem, values) = problem.rescale(values)
