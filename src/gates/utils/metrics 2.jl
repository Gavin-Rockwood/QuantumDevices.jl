"""
    gate_infidelity(model, gate, target; states=nothing, output_states=states,
                    frame=:lab, include_phases=true, energy_shift=:auto, dense=false, evolution_kwargs=(;))

Compare a target unitary with the evolved gate matrix. With `states=nothing`,
evolve the full retained-space identity; otherwise evolve only the supplied
orthonormal ket columns. `output_states` defaults to the input states.
`frame=:interaction` removes idle evolution before comparison.
`include_phases=false` compares probabilities using [`unitary_infidelity`](@ref).
Projected columns retain leakage. Automatic energy centering removes a constant
idle offset during integration and restores global phase. `energy_shift=0`
disables it; explicit offsets use frequency units.
`dense=true` materializes all Hamiltonian matrices during construction.
`evolution_kwargs` are native solver options.
This is the default calibration objective; it does not model dissipation.
"""
function gate_infidelity(model::DeviceModel, gate::DeviceGate, target;
    states=nothing, output_states=states, frame=:lab, include_phases::Bool=true,
    energy_shift=:auto, dense::Bool=false, evolution_kwargs=(;))
    objective = _prepared_gate_objective(model, gate, target;
        states, output_states, frame, include_phases, energy_shift, dense, evolution_kwargs)
    return objective(gate)
end

# This cache belongs to one optimization problem; there is no global model cache.
function _prepared_gate_numerical(model, gate; scalar::Number=1, energy_shift=0, dense::Bool=false)
    _gate_parameters(model, gate) # Validate endpoints and required parameters now.
    dims = [get(model.truncation_dimensions, c, c.dimension) for c in model.components]
    terms = vcat(model.compiled_terms, _projected_terms(gate.hamiltonian, model.components, dims))
    variable = []
    fixed_terms = []
    gate_keys = Set(keys(gate.parameters))
    for term in terms
        destination = isdisjoint(parameter_keys(term[2]), gate_keys) ? fixed_terms : variable
        push!(destination, term)
    end
    fixed = _evaluate_terms(fixed_terms, model.parameters; scalar, dense)
    if !iszero(energy_shift)
        offset = (scalar * energy_shift) * qeye_like(model.H)
        fixed = fixed === nothing ? -offset : fixed - offset
    end
    return candidate -> _evaluate_terms(variable, _gate_parameters(model, candidate); fixed, scalar, dense)
end

function _prepared_gate_objective(model, gate, target;
    states, output_states, frame, include_phases, energy_shift=:auto, dense::Bool=false, evolution_kwargs)
    if states === nothing
        output_states === nothing || throw(ArgumentError("output_states requires input states"))
        initial = qeye_like(model.H)
        output = Matrix(initial.data)
    else
        initial, output = _gate_bases(states, output_states)
        initial.dimensions.to == model.H.dimensions.to ||
            throw(DimensionMismatch("Input basis and model have different dimensions"))
    end
    target_matrix = _validated_unitary_target(target)
    size(target_matrix) == (size(output, 2), size(initial, 2)) ||
        throw(DimensionMismatch("Target must act on the selected gate basis"))
    projection = _prepare_gate_projection(output, frame, model.H)
    shift = _energy_shift(model.H, initial, energy_shift)
    evaluate = _prepared_gate_numerical(model, gate; scalar=2pi, energy_shift=shift, dense)
    return function (candidate)
        options = merge((; tstops=pulse_tstops(candidate)), (; evolution_kwargs...))
        columns = _evolve_gate_columns(evaluate(candidate), initial, candidate.duration; phase_shift=2pi*shift, options...)
        actual = _project_gate_columns(columns, projection, candidate.duration)
        return 1 - _validated_matrix_fidelity(target_matrix, actual, include_phases)
    end
end
