"""Return a physical full-space propagator. Selected-column results have no unitary."""
function unitary(result::GateEvolution; step=length(result.times))
    isnothing(result.inputs) || throw(ArgumentError("Selected-column evolution does not contain a full unitary"))
    return result.propagators[step]
end

"""Unwrapped reference idle phases `2π E t`; idle evolution is `exp(-im*phase)`."""
function idle_phases(result::GateEvolution; step=length(result.times))
    phases = 2pi*result.times[step].*result.reference.energies
    return TrackedData(phases,result.reference.tracking)
end

"""
    gate_matrix(result; states=nothing, output_states=nothing, step=...)

Project physical evolution onto ordered orthonormal bases, preserving leakage.
For a full propagator, `states` selects input kets/labels (energy-ordered reference eigenstates by
default), and output states default to inputs. Selected-column propagation keeps
its recorded input order, rejects `states`, and defaults outputs to stored inputs.
Distinct output bases may have a different number of states.
All physical dynamical phases are retained; projection applies no reference-Hamiltonian
inverse or phase compensation.
"""
function gate_matrix(result::GateEvolution; states=nothing, output_states=nothing,
    step=length(result.times))
    if isnothing(result.inputs)
        inputs = isnothing(states) ? deepcopy(_reference_states(result.reference)) : _gate_kets(result.reference,states)
        columns = result.propagators[step].data * _gate_columns(inputs).data
    else
        isnothing(states) || throw(ArgumentError("Selected-column evolution uses its recorded inputs; use output_states to select outputs"))
        inputs = result.inputs
        columns = result.propagators[step].data
    end
    outputs = isnothing(output_states) ? inputs : _gate_kets(result.reference,output_states)
    return QuantumObject(_gate_columns(outputs).data'*columns;
        dims=QuantumToolbox.Dimensions(QuantumToolbox.Space(length(outputs)),QuantumToolbox.Space(length(inputs))))
end
