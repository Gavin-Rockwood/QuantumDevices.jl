"""
    Measurement(effects::AbstractVector; labels=eachindex(effects), acts_on=nothing, atol=1e-8)
    Measurement(observable::AbstractQuantumObject; acts_on=nothing, atol=1e-8)

An ordered, complete POVM. Effects must be Hermitian positive operators summing
to identity. An observable is converted to projectors onto its distinct
eigenspaces, ordered and labeled by eigenvalue (degeneracy tolerance `atol`).
`acts_on` selects ascending, unique, one-based tensor subsystem indices;
otherwise effects act on the full retained state. Leakage requires an explicit
outcome. Outcome labels describe vector order, not array indices.
"""
struct Measurement{E,L,A,T}
    effects::E
    labels::L
    acts_on::A
    atol::T
    function Measurement(effects::AbstractVector; labels=collect(eachindex(effects)), acts_on=nothing, atol=1e-8)
        isfinite(atol) && atol > 0 || throw(ArgumentError("Measurement tolerance must be finite and positive"))
        isempty(effects) && throw(ArgumentError("A measurement needs outcome effects"))
        effects = copy(effects)
        labels = collect(labels)
        length(labels) == length(effects) && allunique(labels) ||
            throw(ArgumentError("Measurement labels must be unique and match the effects"))
        all(e -> e isa AbstractQuantumObject && isoper(e), effects) ||
            throw(ArgumentError("Measurement effects must be quantum operators"))
        first_effect = first(effects)
        for effect in effects
            effect.dimensions == first_effect.dimensions || throw(DimensionMismatch("Effect dimensions differ"))
            _tomography_positive(effect, atol, "Measurement effect")
        end
        norm(Matrix(sum(effects).data - qeye_like(first_effect).data), Inf) <= atol ||
            throw(ArgumentError("Measurement effects must sum to identity; include leakage outcomes"))
        if acts_on !== nothing
            acts_on = collect(acts_on)
            !isempty(acts_on) && all(i -> i isa Integer && i > 0, acts_on) &&
                issorted(acts_on) && allunique(acts_on) ||
                throw(ArgumentError("acts_on must contain ascending unique positive subsystem indices"))
            acts_on = Tuple(acts_on)
        end
        new{typeof(effects),typeof(labels),typeof(acts_on),typeof(atol)}(effects, labels, acts_on, atol)
    end
end

function _tomography_positive(operator, atol, name)
    matrix = Matrix(operator.data)
    all(isfinite, matrix) || throw(ArgumentError("$name must be finite"))
    norm(matrix - matrix', Inf) <= atol || throw(ArgumentError("$name must be Hermitian"))
    minimum(eigvals(Hermitian(matrix))) >= -atol || throw(ArgumentError("$name must be positive"))
    return nothing
end

function Measurement(observable::AbstractQuantumObject; acts_on=nothing, atol=1e-8)
    isoper(observable) || throw(ArgumentError("Observable must be an operator"))
    matrix = Matrix(observable.data)
    isfinite(atol) && atol > 0 || throw(ArgumentError("Measurement tolerance must be finite and positive"))
    all(isfinite, matrix) && norm(matrix - matrix', Inf) <= atol ||
        throw(ArgumentError("Observable must be finite and Hermitian"))
    values, states = eigenstates(observable)
    order = sortperm(real.(values))
    labels = Float64[]
    effects = typeof(ket2dm(states[first(order)]))[]
    for i in order
        value = real(values[i])
        if !isempty(labels) && abs(value - last(labels)) <= atol
            effects[end] = effects[end] + ket2dm(states[i])
        else
            push!(labels, value)
            push!(effects, ket2dm(states[i]))
        end
    end
    return Measurement(effects; labels, acts_on, atol)
end
