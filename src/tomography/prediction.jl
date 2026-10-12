_tomography_resolve(value, model) = applicable(value, model) ? value(model) : value

function _tomography_scalar_predictions(model, data)
    return map(data) do datum
        prediction = datum.prediction(model)
        prediction isa Real && !(prediction isa Bool) && isfinite(prediction) ||
            throw(ArgumentError("Model prediction must be a finite real scalar"))
        prediction
    end
end

function _tomography_state(state, model; atol=1e-8)
    state isa AbstractQuantumObject && (isket(state) || isoper(state)) ||
        throw(ArgumentError("Preparation must resolve to a ket or density matrix"))
    density = isket(state) ? ket2dm(state) : state
    size(state, 1) == size(model.H, 1) && density.dimensions == model.H.dimensions ||
        throw(DimensionMismatch("Preparation must match the retained model space"))
    all(isfinite, state.data) || throw(ArgumentError("Preparation must be finite"))
    if isket(state)
        isapprox(norm(state), 1; atol, rtol=0) || throw(ArgumentError("Preparation ket must be normalized"))
    else
        _tomography_positive(state, atol, "Preparation density matrix")
        abs(tr(state) - 1) <= atol || throw(ArgumentError("Preparation density matrix must have unit trace"))
    end
    return state
end

function _tomography_probabilities(state, measurement::Measurement; atol=1e-6)
    reduced = measurement.acts_on === nothing ? state : ptrace(isket(state) ? ket2dm(state) : state, measurement.acts_on)
    dimensions = isket(reduced) ? ket2dm(reduced).dimensions : reduced.dimensions
    dimensions == first(measurement.effects).dimensions ||
        throw(DimensionMismatch("Measurement effects must match the measured subsystems"))
    probabilities = [real(expect(effect, reduced)) for effect in measurement.effects]
    all(isfinite, probabilities) && all(p -> -atol <= p <= 1 + atol, probabilities) &&
        abs(sum(probabilities) - 1) <= atol ||
        throw(ArgumentError("Predicted probabilities violate positivity or normalization"))
    probabilities = clamp.(probabilities, 0, 1)
    return probabilities ./ sum(probabilities)
end

function _tomography_propagators(model, sequences; evolution_kwargs, threaded::Bool=false)
    names = unique([name for sequence in sequences for name in sequence])
    for name in names
        haskey(model.gates, name) || throw(ArgumentError("Unknown gate: $name"))
    end
    # Workers write distinct vector slots; construct the dictionary after joining.
    matrices = Vector{Any}(undef, length(names))
    build = function(i)
        gate = model.gates[names[i]]
        matrices[i] = gate.duration == 0 ? qeye_like(model.H) :
            get_unitary(model, gate; evolution_kwargs...)
    end
    if threaded
        Threads.@threads for i in eachindex(names)
            build(i)
        end
    else
        for i in eachindex(names)
            build(i)
        end
    end
    return Dict{Any,Any}(name => matrices[i] for (i, name) in enumerate(names))
end

function _tomography_predict(model, sequence, preparation, measurement, propagators,
    preparations=IdDict(), measurements=IdDict())
    state = get!(preparations, preparation) do
        _tomography_state(_tomography_resolve(preparation, model), model)
    end
    measurement = get!(measurements, measurement) do
        resolved = _tomography_resolve(measurement, model)
        resolved isa Measurement || throw(ArgumentError("Measurement must resolve to Measurement"))
        resolved
    end
    for name in sequence
        U = propagators[name]
        state = isket(state) ? U * state : U * state * U'
    end
    return _tomography_probabilities(state, measurement)
end

"""
    predict_sequence(model, sequence, preparation, measurement;
        evolution_kwargs=(; abstol=1e-8, reltol=1e-8), threaded::Bool=false)

Predict ordered outcome probabilities without observed data. Preparations and
measurements may be fixed objects or functions of the candidate model. Gates use
their local time origins, full retained-space unitaries, and existing frequency
units. Empty sequences measure the preparation; zero-duration gates are identity.
Native solver options pass through `evolution_kwargs`. Set `threaded=true` to
build distinct gate propagators in parallel on Julia threads. Shared controls and
solver callbacks must be safe for concurrent calls.
"""
function predict_sequence(model::DeviceModel, sequence, preparation, measurement;
    evolution_kwargs=(; abstol=1e-8, reltol=1e-8), threaded::Bool=false)
    sequence isa AbstractString && throw(ArgumentError("Pass a collection of gate names"))
    sequence = collect(sequence)
    sequence isa AbstractVector || throw(ArgumentError("Gate sequence must be one-dimensional"))
    propagators = _tomography_propagators(model, (sequence,); evolution_kwargs, threaded)
    return _tomography_predict(model, sequence, preparation, measurement, propagators)
end

"""
    predict_sequences(model, data; evolution_kwargs=(; abstol=1e-8, reltol=1e-8), threaded::Bool=false)

Predict probability vectors in dataset order, ignoring observed distributions.
Each required gate propagator is computed once per call and reused across records.
Repeated preparation and measurement objects are resolved once per call by identity.
No caches survive between candidate models. Inputs are not modified.
`threaded=true` parallelizes propagator construction only; shared controls and
solver callbacks must be safe for concurrent calls.
"""
function predict_sequences(model::DeviceModel, data; evolution_kwargs=(; abstol=1e-8, reltol=1e-8), threaded::Bool=false)
    data = collect(data)
    all(d -> d isa SequenceData, data) || throw(ArgumentError("Expected SequenceData records"))
    propagators = _tomography_propagators(model, (d.sequence for d in data); evolution_kwargs, threaded)
    preparations, measurements = IdDict(), IdDict()
    return [_tomography_predict(model, d.sequence, d.preparation, d.measurement,
        propagators, preparations, measurements) for d in data]
end
