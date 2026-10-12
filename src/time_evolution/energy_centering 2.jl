# Energy offsets are real constants, not a rotating-wave approximation.
# Helpers use the units of the supplied Hamiltonian (angular units for sesolve).
function _energy_shift(H::QuantumObject, initial::QuantumObject, option)
    if option === :auto
        columns = isket(initial) ? reshape(initial.data, :, 1) : initial.data
        action = H.data * columns
        energy = 0.0
        count = 0
        for j in axes(columns, 2)
            state = view(columns, :, j)
            weight = real(dot(state, state))
            iszero(weight) && continue
            count += 1
            energy += real(dot(state, view(action, :, j))) / weight
        end
        return iszero(count) ? 0.0 : energy / count
    end
    option isa Real && isfinite(option) || throw(ArgumentError(
        "energy_shift must be :auto or a finite real energy offset (0 disables centering)"))
    return option
end

function _restore_global_phase!(solution::QuantumToolbox.TimeEvolutionSol, shift, start)
    iszero(shift) && return solution
    for (state, time) in zip(solution.states, solution.times_states)
        rmul!(state.data, cis(-shift * (time - start)))
    end
    return solution
end

function _sesolve_centered(H, initial, times; energy_shift=:auto, kwargs...)
    isempty(times) && throw(ArgumentError("Evolution times cannot be empty"))
    hamiltonian = H isa Tuple ? QobjEvo(H) : H
    params = get(kwargs, :params, SciMLBase.NullParameters())
    snapshot = hamiltonian isa QobjEvo ? hamiltonian(params, first(times)) : hamiltonian
    shift = _energy_shift(snapshot, initial, energy_shift)
    centered = iszero(shift) ? hamiltonian : hamiltonian - shift * qeye_like(snapshot)
    solution = QuantumToolbox.sesolve(centered, initial, times; kwargs...)
    return _restore_global_phase!(solution, shift, first(times))
end
