# Capture the fixed Hamiltonian and input/output basis once. Each invocation owns
# its runtime pulses and QobjEvo workspace; no reference or compiler work is repeated.
function _calibration_evaluator(prepared,reference,selections;
    paths=nothing,remove_vacuum::Bool=false,energy_shift=0,herm_check=1e-9,kwargs...)
    inputs = _gate_kets(reference,selections)
    initial = _gate_columns(inputs)
    projection = initial.data'
    vacuum_energy = remove_vacuum ? reference_vacuum_energy(prepared.model) : 0
    shift = if energy_shift === :auto
        coefficients = _reference_eigenvectors(reference).data' * initial.data
        sum(reference.energies[i]*abs2(coefficients[i,j])
            for i in axes(coefficients,1),j in axes(coefficients,2))/size(coefficients,2)-vacuum_energy
    else
        energy_shift isa Real && isfinite(energy_shift) ||
            throw(ArgumentError("energy_shift must be :auto or a finite real GHz offset"))
        energy_shift
    end
    hamiltonian = 2pi*_compile_gate(prepared.model,prepared.recipe,prepared.terms,prepared.parameters;
        remove_vacuum,herm_check,runtime=isnothing(paths),free_parameters=isnothing(paths) ? () : paths)
    iszero(shift) || (hamiltonian = hamiltonian-2pi*shift*_reference_identity(reference))
    dimensions = QuantumToolbox.Dimensions(QuantumToolbox.Space(length(inputs)),QuantumToolbox.Space(length(inputs)))
    prepare=isnothing(paths) ? nothing : _compile_gate_preparation(prepared.recipe,prepared.parameters,paths)
    return function(candidate)
        prepared_runtime = isnothing(paths) ? nothing : prepare(candidate)
        recipe = isnothing(paths) ? candidate : prepared.recipe
        duration = isnothing(paths) ? recipe.duration : prepared_runtime.duration
        isfinite(duration) && duration>0 || throw(ArgumentError("Evolution duration must be finite and positive"))
        pulses = isnothing(paths) ? Dict(c.target=>let
            resolved_recipe=_resolve_pulse(c.pulse.recipe,prepared.parameters)
            waveform=_compile_waveform(resolved_recipe)
            t -> waveform(nothing,t,t)
        end for c in recipe.controls) : prepared_runtime.runtime
        options = _gate_solver_options(recipe,duration;runtime_stops=isnothing(paths) ? nothing : prepared_runtime.stops,kwargs...)
        solution = sesolve(deepcopy(hamiltonian),copy(initial),[zero(duration),duration];params=pulses,options...)
        SciMLBase.successful_retcode(solution.retcode) || error("Gate evolution failed with return code $(solution.retcode)")
        isempty(solution.states) && error("Gate evolution saved no propagators")
        last(solution.times_states)==duration || error("Gate evolution did not save the final propagator")
        columns = last(solution.states).data
        matrix = projection*columns
        iszero(shift) || rmul!(matrix,cis(-2pi*shift*duration))
        return QuantumObject(matrix;dims=dimensions)
    end
end
