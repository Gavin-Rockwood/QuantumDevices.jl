"""
Physical lab-frame evolution in the reference model's construction coordinates. `inputs === nothing` denotes full identity propagation; otherwise inputs
are ordered orthonormal kets. `times` are elapsed times and always include the
final gate duration. Propagators retain leakage and all relative physical phases,
with an explicitly selected vacuum energy origin recorded in solver metadata.

Reference arrays and provenance are captured independently of the prepared model.
As for GateRecipe, arbitrary mutable state hidden inside user closures is outside
snapshot guarantees. Arrays and QuantumObjects retain ordinary mutable behavior.
"""
struct GateEvolution{G,R,I,T,P,S}
    recipe::G
    reference::R
    inputs::I
    times::T
    propagators::P
    solver::S
end

function _gate_reference(model)
    source = model.eigensystem.tracking
    n = length(model.spectrum)
    # Capture actual eigenstates for label inputs, reference-energy diagnostics,
    # and centering expectations in the gate Hamiltonian's construction coordinates.
    tracking = TrackingResult(source.keys, (1,), 1, deepcopy(getfield(source,:_states)), (;),
        copy(getfield(source,:_indices)), deepcopy(getfield(source,:_accuracy)),
        deepcopy(getfield(source,:_diagnostics)); tracked=source.tracked)
    return (; energies=collect(model.spectrum), tracking,
        components=deepcopy(model.eigensystem.components), recipe=deepcopy(model.recipe))
end

function _gate_kets(reference, selections)
    selections isa AbstractVector || throw(ArgumentError("states must be an ordered vector of labels or kets"))
    isempty(selections) && throw(ArgumentError("Pass at least one gate basis state"))
    n = length(reference.energies)
    kets = map(selections) do s
        s isa Tuple ? deepcopy(reference.tracking.states(s)) :
        s isa QuantumObject && isket(s) ? deepcopy(s) :
        throw(ArgumentError("states must contain occupation tuples or QuantumObject kets"))
    end
    all(s -> length(s.data)==n,kets) || throw(DimensionMismatch("Gate kets must match the retained dimension"))
    all(s -> s.dimensions==first(_reference_states(reference)).dimensions,kets) || throw(DimensionMismatch("Gate kets must share quantum dimensions"))
    matrix = hcat((s.data for s in kets)...)
    isapprox(matrix'*matrix, Matrix{ComplexF64}(I,length(kets),length(kets));atol=1e-8,rtol=1e-8) ||
        throw(ArgumentError("Gate basis states must be orthonormal"))
    return kets
end

_gate_columns(kets) = QuantumObject(hcat((s.data for s in kets)...);
    dims=QuantumToolbox.Dimensions(first(kets).dimensions.to, QuantumToolbox.Space(length(kets))))

_reference_states(reference) = getfield(reference.tracking, :_states)
_reference_eigenvectors(reference) = projector_from_states(_reference_states(reference))
function _reference_identity(reference)
    n = length(reference.energies)
    dimensions = Tuple(first(_reference_states(reference)).dims[1])
    return qeye(n; dims=dimensions)
end
