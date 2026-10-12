"""
    state_amplitudes(references::AbstractDict, history::AbstractVector)

Return a dictionary with the same keys as `references`, mapping each reference
state to its complex overlaps `dot(reference, state)` throughout `history`.
The first argument of the overlap is conjugated. Neither references nor saved
states are normalized, and reference states need not be orthogonal or complete.
Use `Dict(k => abs2.(v) for (k, v) in amplitudes)` to obtain squared overlaps;
these are populations when the input states are normalized.

States may be numerical vectors or QuantumToolbox kets. Their dimensions must
match; operators and density matrices do not have state amplitudes. An empty
history produces an empty vector for each reference.
"""
function state_amplitudes(references::AbstractDict, history::AbstractVector)
    vectors = map(_amplitude_vector, history)
    return Dict(key => begin
        reference_vector = _amplitude_vector(reference)
        for (i, state) in enumerate(history)
            length(reference_vector) == length(vectors[i]) || throw(DimensionMismatch(
                "Reference $(repr(key)) and state $i have different lengths"))
            if reference isa AbstractQuantumObject && state isa AbstractQuantumObject
                reference.dimensions == state.dimensions || throw(DimensionMismatch(
                    "Reference $(repr(key)) and state $i have different quantum dimensions"))
            end
        end
        isempty(vectors) ? typeof(complex(zero(eltype(reference_vector))))[] :
            [complex(dot(reference_vector, vector)) for vector in vectors]
    end for (key, reference) in references)
end

function _amplitude_vector(state)
    if state isa AbstractQuantumObject
        isket(state) || throw(ArgumentError("State amplitudes require kets, not operators or density matrices"))
        return state.data
    elseif state isa AbstractVector{<:Number}
        return state
    end
    throw(ArgumentError("Expected a numerical state vector or QuantumToolbox ket, received $(typeof(state))"))
end
