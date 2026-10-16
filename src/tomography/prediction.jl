_tomography_resolve(value, model) = applicable(value, model) ? value(model) : value

_tomography_real_scalar(value) = value isa Real && !(value isa Bool) && isfinite(value)

function _tomography_named_values(values, names, description)
    values isa NamedTuple || throw(ArgumentError("$description must be a NamedTuple"))
    missing = setdiff(names, keys(values))
    extra = setdiff(keys(values), names)
    isempty(missing) && isempty(extra) ||
        throw(ArgumentError("$description names differ: missing $(Tuple(missing)), extra $(Tuple(extra))"))
    return map(names) do name
        value = values[name]
        _tomography_real_scalar(value) ||
            throw(ArgumentError("$description for $name must be a finite real scalar"))
        value
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

function _tomography_probabilities(state, measurement::Measurement; atol=1e-6, leakage_outcome=nothing, normalize_probabilities::Bool=false)
    reduced = measurement.acts_on === nothing ? state : ptrace(isket(state) ? ket2dm(state) : state, measurement.acts_on)
    dimensions = isket(reduced) ? ket2dm(reduced).dimensions : reduced.dimensions
    dimensions == first(measurement.effects).dimensions ||
        throw(DimensionMismatch("Measurement effects must match the measured subsystems"))
    probabilities = [real(expect(effect, reduced)) for effect in measurement.effects]
    mass = sum(probabilities)
    if leakage_outcome !== nothing
        probabilities[leakage_outcome] += max(zero(mass), 1 - mass)
    elseif !normalize_probabilities && mass < 1 - atol
        throw(ArgumentError("Projected evolution lost population; supply a complementary leakage outcome or enlarge subspace"))
    end
    all(isfinite, probabilities) && all(p -> p >= -atol, probabilities) &&
        sum(probabilities) > 0 &&
        (normalize_probabilities || (all(p -> p <= 1 + atol, probabilities) &&
            abs(sum(probabilities) - 1) <= atol)) ||
        throw(ArgumentError("Predicted probabilities violate positivity or normalization: total=$(sum(probabilities)), range=$(extrema(probabilities)), tolerance=$atol"))
    probabilities = max.(probabilities, 0)
    return probabilities ./ sum(probabilities)
end

function _tomography_propagators(model, sequences; evolution_kwargs, threaded::Bool=false, basis=nothing, gates=model.gates, gate_timings=nothing)
    names = unique([name for sequence in sequences for name in sequence])
    for name in names
        haskey(gates, name) || throw(ArgumentError("Unknown gate: $name"))
    end
    # Workers write distinct vector slots; construct the dictionary after joining.
    matrices = Vector{Any}(undef, length(names))
    elapsed = gate_timings === nothing ? nothing : zeros(length(names))
    build = function(i)
        started = time_ns()
        gate = gates[names[i]]
        matrices[i] = if basis === nothing
            gate.duration == 0 ? qeye_like(model.H) : get_unitary(model, gate; evolution_kwargs...)
        else
            gate.duration == 0 ? qeye(length(basis)) :
                numerical(model, gate, basis; frame=:lab, evolution_kwargs...)
        end
        elapsed === nothing || (elapsed[i] = (time_ns() - started) / 1e9)
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
    gate_timings === nothing || append!(gate_timings, [names[i] => elapsed[i] for i in eachindex(names)])
    return Dict{Any,Any}(name => matrices[i] for (i, name) in enumerate(names))
end

function _tomography_predict(model, sequence, preparation, measurement, propagators,
    preparations=IdDict(), measurements=IdDict(); basis=nothing, leakage_outcome=nothing, normalize_probabilities::Bool=false)
    state = get!(preparations, preparation) do
        _tomography_state(_tomography_resolve(preparation, model), model)
    end
    measurement = get!(measurements, measurement) do
        resolved = _tomography_resolve(measurement, model)
        resolved isa Measurement || throw(ArgumentError("Measurement must resolve to Measurement"))
        resolved
    end
    columns = basis === nothing ? nothing : _orthonormal_columns(basis)
    if columns !== nothing
        state = isket(state) ? columns' * state : columns' * state * columns
        mass = isket(state) ? real(dot(state,state)) : real(tr(state))
        isapprox(mass,1; atol=1e-8,rtol=0) ||
            throw(ArgumentError("Preparation must lie in the supplied subspace"))
        if leakage_outcome !== nothing
            leakage_outcome isa Integer && !(leakage_outcome isa Bool) &&
                1 <= leakage_outcome <= length(measurement.effects) ||
                throw(ArgumentError("leakage_outcome must index a measurement outcome"))
            for ket in basis
                _tomography_probabilities(ket,measurement)[leakage_outcome] <= 1e-8 ||
                    throw(ArgumentError("The leakage effect must annihilate every subspace state"))
            end
            outside = qeye_like(model.H) - columns * columns'
            reduced_outside = measurement.acts_on === nothing ? outside : ptrace(outside,measurement.acts_on)
            outside_rank = size(model.H,1) - length(basis)
            isapprox(real(tr(measurement.effects[leakage_outcome]*reduced_outside)),outside_rank;
                atol=1e-8,rtol=0) || throw(ArgumentError("The leakage effect must include the entire complement of subspace"))
        end
    elseif leakage_outcome !== nothing
        throw(ArgumentError("leakage_outcome requires subspace"))
    end
    for name in sequence
        U = propagators[name]
        state = isket(state) ? U * state : U * state * U'
        if columns === nothing && !normalize_probabilities
            # Full closed-system evolution preserves mass. Correct small solver
            # drift per gate before it accumulates across repeated sequences.
            # Projected evolution loses real population and must not be normalized.
            mass = isket(state) ? real(dot(state, state)) : real(tr(state))
            isfinite(mass) && abs(mass - 1) <= 1e-6 ||
                throw(ArgumentError("Gate $name changed state normalization to $mass; tighten evolution_kwargs abstol/reltol"))
            state = state / (isket(state) ? sqrt(mass) : mass)
        end
    end
    if columns !== nothing
        state = isket(state) ? columns * state : columns * state * columns'
    end
    # A requested normalization must not hide discarded projected population.
    normalize = normalize_probabilities && (columns === nothing || leakage_outcome !== nothing)
    return _tomography_probabilities(state, measurement; leakage_outcome, normalize_probabilities=normalize)
end

"""
    predict_sequence(model, sequence; preparation, measurement, kwargs...)
    predict_sequence(model, sequence, preparation, measurement;
        evolution_kwargs=(; abstol=1e-8, reltol=1e-8), threaded::Bool=false, subspace=nothing, leakage_outcome=nothing, normalize_probabilities::Bool=false, gates=model.gates)

Predict ordered outcome probabilities without observed data. Preparations and
measurements may be fixed objects or functions of the candidate model. Gates use
their local time origins, full retained-space unitaries, and existing frequency
units. Empty sequences measure the preparation; zero-duration gates are identity.
`subspace` may be an ordered orthonormal ket vector or `model -> vector`.
Only those columns evolve in the full Hamiltonian; each gate is projected to the
basis before multiplying sequences. This approximation discards leakage return.
Preparations must lie in the basis. To preserve discarded probability, supply
`leakage_outcome` indexing the complementary measurement effect that annihilates
the basis. Lost population is an absorbing outcome, never renormalized away.
`normalize_probabilities=true` normalizes finite nonnegative outcome weights at
measurement, after assigning discarded population to the complementary leakage
outcome. It does not normalize projected states or permit unaccounted population
loss when no leakage outcome is supplied. Negative weights beyond rounding
tolerance still raise an error. The default `false` retains strict checks.
With `normalize_probabilities=false`, full-space evolution corrects normalization
drift of at most `1e-6` after each
gate; larger drift is rejected. Projected population loss is preserved.
`gates` explicitly selects the gate dictionary (default: `model.gates`).
Native solver options pass through `evolution_kwargs`. Set `threaded=true` to
build distinct gate propagators in parallel on Julia threads. Shared controls and
solver callbacks must be safe for concurrent calls.
"""
function predict_sequence(model::DeviceModel, sequence, preparation, measurement;
    evolution_kwargs=(; abstol=1e-8, reltol=1e-8), threaded::Bool=false, subspace=nothing, leakage_outcome=nothing, normalize_probabilities::Bool=false, gates=model.gates)
    sequence isa AbstractString && throw(ArgumentError("Pass a collection of gate names"))
    sequence = collect(sequence)
    sequence isa AbstractVector || throw(ArgumentError("Gate sequence must be one-dimensional"))
    basis = subspace === nothing ? nothing : _tomography_resolve(subspace,model)
    if basis !== nothing
        columns = _orthonormal_columns(basis)
        columns.dimensions.to == model.H.dimensions.to ||
            throw(DimensionMismatch("Subspace must match retained model dimensions"))
    end
    propagators = _tomography_propagators(model, (sequence,); evolution_kwargs, threaded, basis, gates)
    return _tomography_predict(model, sequence, preparation, measurement, propagators; basis, leakage_outcome, normalize_probabilities)
end

predict_sequence(model::DeviceModel, sequence; preparation, measurement, kwargs...) =
    predict_sequence(model,sequence,preparation,measurement;kwargs...)

"""
    predict_sequences(model, sequences; preparation, measurement, kwargs...)
    predict_sequences(model, data; evolution_kwargs=(; abstol=1e-8, reltol=1e-8), threaded::Bool=false, subspace=nothing, leakage_outcome=nothing, normalize_probabilities::Bool=false, gates=model.gates)

Predict directly from a list of gate sequences using shared preparation and
measurement objects or callbacks; no observed data is required. Alternatively,
pass `SequenceData` records; their observed distributions are ignored.
Return probability vectors in input sequence order.
Each required gate propagator is computed once per call and reused across records.
Repeated preparation and measurement objects are resolved once per call by identity.
`subspace` and `leakage_outcome` follow [`predict_sequence`](@ref): projected
matrices are computed once per required gate; discarded population cannot return.
No caches survive between candidate models. Inputs are not modified.
Optional `diagnostics=report -> ...` receives seconds for `setup_seconds`,
`propagation_seconds`, `sequence_seconds`, `total_seconds`, and `gate_seconds`
(gate-name => seconds pairs). Threaded gate times overlap. Callback execution
is excluded from the total; callback errors propagate.
`threaded=true` parallelizes propagator construction only; shared controls and
solver callbacks must be safe for concurrent calls.
"""
function predict_sequences(model::DeviceModel, data; preparation=nothing, measurement=nothing, evolution_kwargs=(; abstol=1e-8, reltol=1e-8), threaded::Bool=false, subspace=nothing, leakage_outcome=nothing, normalize_probabilities::Bool=false, gates=model.gates, diagnostics=nothing)
    started = time_ns()
    data = collect(data)
    if preparation === nothing && measurement === nothing
        all(d -> d isa SequenceData, data) ||
            throw(ArgumentError("Plain sequences require preparation and measurement"))
    else
        preparation !== nothing && measurement !== nothing ||
            throw(ArgumentError("Supply both preparation and measurement"))
        any(d -> d isa SequenceData, data) && throw(ArgumentError("Do not mix SequenceData with shared preparation/measurement"))
        data = [(;sequence=s,preparation,measurement) for s in _tomography_sequences(data)]
    end
    basis = subspace === nothing ? nothing : _tomography_resolve(subspace,model)
    if basis !== nothing
        columns = _orthonormal_columns(basis)
        columns.dimensions.to == model.H.dimensions.to ||
            throw(DimensionMismatch("Subspace must match retained model dimensions"))
    end
    setup_seconds = (time_ns() - started) / 1e9
    gate_timings = diagnostics === nothing ? nothing : Pair{Any,Float64}[]
    propagation_started = time_ns()
    propagators = _tomography_propagators(model, (d.sequence for d in data); evolution_kwargs, threaded, basis, gates, gate_timings)
    propagation_seconds = (time_ns() - propagation_started) / 1e9
    prediction_started = time_ns()
    preparations, measurements = IdDict(), IdDict()
    result = [_tomography_predict(model, d.sequence, d.preparation, d.measurement,
        propagators, preparations, measurements; basis, leakage_outcome, normalize_probabilities) for d in data]
    diagnostics === nothing || diagnostics((; setup_seconds, propagation_seconds,
        sequence_seconds=(time_ns()-prediction_started)/1e9, gate_seconds=gate_timings,
        total_seconds=(time_ns()-started)/1e9))
    return result
end
