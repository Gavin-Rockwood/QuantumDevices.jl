"""
    ModelTomographyProblem(model, data, paths; loss=:auto,
        probability_floor=1e-12, evolution_kwargs=(; abstol=1e-8, reltol=1e-8), threaded::Bool=false, kwargs...)

Fit finite real model/gate scalars selected by writable slash paths. Accepts
either [`SequenceData`](@ref) records or [`ModelDatum`](@ref) scalar observations;
the two types cannot be mixed in one fit. Exposes a
standard SciML optimization `.problem`, immutable `.rebuild(values)`, and `.paths`.
Solve using `solve(setup, algorithm; kwargs...)`; no optimizer is selected here.
Candidate values are applied together: each changed component is reconstructed
once and dressed states are retracked once when physical parameters change.
Gate-only candidates reuse the idle model. Every candidate starts from the input
model, independently of earlier optimizer evaluations. Paths must resolve to
distinct writable destinations; aliases and overlapping paths are rejected.
`loss=:auto` selects `:likelihood` when every record has shots and `:least_squares`
when none does. Mixed shot-qualified sequence datasets require explicit
`:least_squares`. Scalar datasets always use least squares, with residuals
`(prediction(model) - value) / uncertainty`. They reconstruct each candidate
once and evaluate all predictions without gate evolution.
Likelihood is summed multinomial negative log likelihood without constants;
probabilities with shots use `-N*sum(f .* log(p))`, without rounding counts.
Only logarithms use `probability_floor`; zero-count terms contribute zero.
Least squares sums squared outcome-probability residuals without observable-
variance weighting. Additional keywords, including bounds, pass to SciML.
`threaded=true` builds distinct gate propagators in parallel on Julia threads;
it and `evolution_kwargs` have no effect on scalar fits.
Shared controls and solver callbacks must be safe for concurrent calls.
Automatic differentiation is not guaranteed through model reconstruction.
"""
struct ModelTomographyProblem{P,F,R,S}
    problem::P
    rebuild::F
    residuals::R
    paths::S
end

function _tomography_parameter_setup(model, paths)
    paths isa AbstractString && throw(ArgumentError("Pass a collection of parameter paths"))
    paths = collect(paths)
    !isempty(paths) && all(p -> p isa AbstractString, paths) && allunique(paths) ||
        throw(ArgumentError("Parameter paths must be a nonempty unique collection of strings"))
    update_tree = _prepare_model_updates(model, paths)
    initial = Float64[]
    for path in paths
        value = getpath(model, path)
        value isa Real && isfinite(value) || throw(ArgumentError("Parameter $path must be finite and real"))
        push!(initial, value)
    end
    rebuild = function(values)
        length(values) == length(paths) || throw(DimensionMismatch("Candidate length does not match selected paths"))
        all(v -> v isa Real && isfinite(v), values) || throw(ArgumentError("Candidate values must be finite real scalars"))
        return _apply_model_updates(model, update_tree, values)
    end
    return paths, initial, rebuild
end

function ModelTomographyProblem(model::DeviceModel, data, paths; loss=:auto,
    probability_floor=1e-12, evolution_kwargs=(; abstol=1e-8, reltol=1e-8), threaded::Bool=false, kwargs...)
    paths, initial, rebuild = _tomography_parameter_setup(model, paths)
    data = collect(data)
    !isempty(data) || throw(ArgumentError("A tomography problem needs observations"))
    if all(d -> d isa ModelDatum, data)
        loss in (:auto, :least_squares) || throw(ArgumentError("ModelDatum observations require least squares"))
        residuals = function(values)
            predictions = _tomography_scalar_predictions(rebuild(values), data)
            return [(prediction - datum.value) / datum.uncertainty
                for (prediction, datum) in zip(predictions, data)]
        end
        residuals(initial)
        objective = (values, _) -> sum(abs2, residuals(values))
        problem = SciMLBase.OptimizationProblem(objective, initial; kwargs...)
        return ModelTomographyProblem(problem, rebuild, residuals, paths)
    end
    all(d -> d isa SequenceData, data) ||
        throw(ArgumentError("Use only ModelDatum or only SequenceData records; fit scalar constraints first, then refine with sequence data"))
    qualified = [d.shots !== nothing for d in data]
    if loss === :auto
        loss = all(qualified) ? :likelihood : !any(qualified) ? :least_squares :
            throw(ArgumentError("Mixed shot-qualified data require loss=:least_squares"))
    end
    loss in (:likelihood, :least_squares) || throw(ArgumentError("Unknown tomography loss: $loss"))
    loss === :likelihood && !all(qualified) && throw(ArgumentError("Likelihood requires shots for every record"))
    isfinite(probability_floor) && 0 < probability_floor < 1 || throw(ArgumentError("Probability floor must lie strictly between zero and one"))
    predictions = values -> predict_sequences(rebuild(values), data; evolution_kwargs, threaded)
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
    problem = SciMLBase.OptimizationProblem(objective, initial; kwargs...)
    return ModelTomographyProblem(problem, rebuild, loss === :least_squares ? residuals : nothing, paths)
end

"""
    reconstructed_model(problem::ModelTomographyProblem, values)

Rebuild a fitted model in selected parameter order without modifying the input.
Inspect the optimizer retcode and fit quality before accepting its values.
"""
reconstructed_model(problem::ModelTomographyProblem, values) = problem.rebuild(values)

"""
    tomography_residuals(problem::ModelTomographyProblem, values)

Return prediction-minus-observation residuals in dataset and outcome order.
For [`ModelDatum`](@ref), each residual is divided by its uncertainty.
Only available for least squares; their squared norm equals the scalar objective.
External residual-based solvers may call this function directly.
"""
function tomography_residuals(problem::ModelTomographyProblem, values)
    problem.residuals === nothing && throw(ArgumentError("Likelihood problems do not have least-squares residuals"))
    return problem.residuals(values)
end
