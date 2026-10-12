"""
    ModelDatum(prediction, value; uncertainty=1)

One scalar model observation. `prediction(model)` must return a finite real
scalar, such as a difference of dressed energies in `model.spectrum`. `value`
must be finite and real; `uncertainty` must also be positive. Observations,
predictions, and uncertainties use the same units. The least-squares residual
is `(prediction(model) - value) / uncertainty`; the default is unweighted.
Errors are treated as independent. No gate evolution is needed for scalar fits.
"""
struct ModelDatum{F,V,U}
    prediction::F
    value::V
    uncertainty::U
    function ModelDatum(prediction, value; uncertainty=1)
        value isa Real && !(value isa Bool) && isfinite(value) ||
            throw(ArgumentError("Model observation must be a finite real scalar"))
        uncertainty isa Real && !(uncertainty isa Bool) && isfinite(uncertainty) && uncertainty > 0 ||
            throw(ArgumentError("Model uncertainty must be a finite positive real scalar"))
        new{typeof(prediction),typeof(value),typeof(uncertainty)}(prediction, value, uncertainty)
    end
end

"""
    SequenceData(sequence, preparation, measurement; counts=nothing,
                 probabilities=nothing, shots=nothing)

Measured outcomes following named gates in execution order. Supply a normalized
ket/density matrix or `model -> state` preparation, and a [`Measurement`](@ref)
or `model -> Measurement`. Exactly one distribution input is required, in the
measurement's outcome order. Counts are nonnegative integers with positive total;
shots are inferred. Probabilities may optionally carry a positive integer shot
total, and are preserved without rounding into counts. Input vectors are copied.
"""
struct SequenceData{S,P,M,F,C,N}
    sequence::S
    preparation::P
    measurement::M
    probabilities::F
    counts::C
    shots::N
    function SequenceData(sequence, preparation, measurement; counts=nothing, probabilities=nothing, shots=nothing)
        sequence isa AbstractString && throw(ArgumentError("Pass a collection of gate names"))
        (counts === nothing) != (probabilities === nothing) ||
            throw(ArgumentError("Supply exactly one of counts or probabilities"))
        shots === nothing || (shots isa Integer && !(shots isa Bool) && shots > 0) ||
            throw(ArgumentError("Shots must be a positive integer"))
        if counts !== nothing
            counts = collect(counts)
            counts isa AbstractVector && !isempty(counts) && all(n -> n isa Integer && !(n isa Bool) && n >= 0, counts) ||
                throw(ArgumentError("Counts must be nonnegative integers"))
            total = sum(big.(counts))
            0 < total <= typemax(Int) || throw(ArgumentError("Count total must be positive and fit Int"))
            shots === nothing || shots == total || throw(ArgumentError("Shots disagree with counts"))
            shots = Int(total)
            probabilities = counts ./ shots
        else
            probabilities = collect(probabilities)
        end
        probabilities isa AbstractVector && !isempty(probabilities) && all(p -> p isa Real && isfinite(p) && 0 <= p <= 1, probabilities) &&
            isapprox(sum(probabilities), 1; atol=1e-8, rtol=0) ||
            throw(ArgumentError("Probabilities must be finite, nonnegative, and sum to one"))
        sequence = collect(sequence)
        sequence isa AbstractVector || throw(ArgumentError("Gate sequence must be one-dimensional"))
        new{typeof(sequence),typeof(preparation),typeof(measurement),typeof(probabilities),typeof(counts),typeof(shots)}(
            sequence, preparation, measurement, probabilities, counts, shots)
    end
end
