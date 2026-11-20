"""
    evolve(prepared::PreparedGate; states=nothing, saveat=nothing, remove_vacuum=false, energy_shift=0, kwargs...)

Evolve the full retained Hamiltonian using `2π * make_hamiltonian(prepared)` and
QuantumToolbox.sesolve. Default propagation is the full identity, saved only at
completion. `states` selects ordered occupation labels (reference eigenstates)
or orthonormal QuantumObject kets. `saveat` requests elapsed-time history; the final
time is appended if necessary. `energy_shift` is a constant offset in GHz (or
`:auto` for the mean reference energy of the input columns). It is subtracted
only during integration; every saved propagator has that centering phase
restored. `remove_vacuum=true` instead leaves the fixed reference vacuum energy
subtracted from the returned evolution, as in `make_hamiltonian`. Automatic
centering is relative to that selected energy origin. The solver metadata records
both choices; neither flag reaches the native ODE solver. Native solver
keywords pass through. Pulse stops are merged with user `tstops`. Unsuccessful
integration throws an error.

Saving must preserve all retained rows/columns: `save_idxs`, `e_ops`, and
`save_end=false` are unsupported. Stored propagators retain all relative physical phases and the explicitly
selected global energy origin. Reference `idle_phases` diagnostics keep the
original reference energies; `result.solver.vacuum_energy` records the removed
global offset.
"""
function evolve(prepared::PreparedGate; states=nothing, saveat=nothing, tstops=Real[], herm_check=1e-9, remove_vacuum::Bool=false, energy_shift=0, kwargs...)
    duration = prepared.recipe.duration
    isfinite(duration) && duration > 0 || throw(ArgumentError("Evolution duration must be finite and positive"))
    for key in (:save_idxs, :e_ops)
        haskey(kwargs,key) && throw(ArgumentError("$key is incompatible with storing physical propagators"))
    end
    get(kwargs,:save_end,true) || throw(ArgumentError("Gate evolution requires save_end=true"))
    reference = _gate_reference(prepared.model)
    inputs = isnothing(states) ? nothing : _gate_kets(reference,states)
    initial = isnothing(inputs) ? _reference_identity(reference) : _gate_columns(inputs)
    vacuum_energy = remove_vacuum ? reference_vacuum_energy(prepared.model) : 0
    shift = if energy_shift === :auto
        coefficients = _reference_eigenvectors(reference).data' * initial.data
        sum(reference.energies[i] * abs2(coefficients[i,j])
            for i in axes(coefficients,1), j in axes(coefficients,2)) / size(coefficients,2) - vacuum_energy
    else
        energy_shift isa Real && isfinite(energy_shift) ||
            throw(ArgumentError("energy_shift must be :auto or a finite real GHz offset"))
        energy_shift
    end
    saves = isnothing(saveat) ? [duration] : saveat isa Real ?
        (isfinite(saveat) && saveat > 0 ? collect(0:saveat:duration) : throw(ArgumentError("saveat spacing must be finite and positive"))) : float.(collect(saveat))
    all(t -> t isa Real && isfinite(t) && 0<=t<=duration,saves) || throw(ArgumentError("saveat times must be finite and within the gate duration"))
    issorted(saves) && allunique(saves) || throw(ArgumentError("saveat times must be increasing and unique"))
    isempty(saves) || last(saves)==duration || push!(saves,duration)
    isempty(saves) && push!(saves,duration)
    stops = collect(tstops)
    all(t -> t isa Real && isfinite(t) && 0<=t<=duration,stops) || throw(ArgumentError("tstops must lie within the gate duration"))
    stops = sort!(unique(vcat(stops,pulse_tstops(prepared))))
    options = merge((;progress_bar=false,save_start=(0 in saves),save_end=true), (;kwargs...), (;saveat=saves,tstops=stops))
    hamiltonian = 2pi*make_hamiltonian(prepared;herm_check,remove_vacuum)
    iszero(shift) || (hamiltonian = hamiltonian - 2pi*shift*_reference_identity(reference))
    solution = sesolve(hamiltonian,initial,[zero(duration),duration];options...)
    SciMLBase.successful_retcode(solution.retcode) || error("Gate evolution failed with return code $(solution.retcode)")
    isempty(solution.states) && error("Gate evolution saved no propagators")
    last(solution.times_states)==duration || error("Gate evolution did not save the final propagator")
    if !iszero(shift)
        for (propagator,time) in zip(solution.states,solution.times_states)
            rmul!(propagator.data,cis(-2pi*shift*time))
        end
    end
    return GateEvolution(deepcopy(prepared.recipe),reference,inputs,collect(solution.times_states),
        solution.states,(;options,remove_vacuum,vacuum_energy,energy_shift=shift,retcode=solution.retcode,alg=solution.alg,
            abstol=solution.abstol,reltol=solution.reltol))
end
