function _floquet_period(T::Real)
    isfinite(T) && T > 0 || throw(ArgumentError("Floquet period must be finite and positive"))
    return T
end

function _floquet_propagator(H, duration, t0; kwargs...)
    initial = qeye_like(H isa QobjEvo ? H(t0) : H)
    isoper(initial) || throw(ArgumentError("Floquet Hamiltonian must be an operator"))
    options = merge((; progress_bar=false), (; kwargs...))
    solution = sesolve(2pi * H, initial, [t0, t0 + duration]; options...)
    SciMLBase.successful_retcode(solution.retcode) ||
        error("Floquet evolution failed with return code $(solution.retcode)")
    isempty(solution.states) && throw(ArgumentError("Floquet evolution requires saved states"))
    return last(solution.states)
end

"""
    get_floquet_basis(H::AbstractQuantumObject, T; t0=0,
                      propagator_kwargs=Dict{Symbol,Any}(), kwargs...)

Diagonalize evolution over `[t0, t0 + T]` for a periodic operator or `QobjEvo`.
Hamiltonians and returned quasienergies use cycles per unit time (GHz for ns).
`basis.modes(t)` returns periodic Floquet kets, so physical solutions are
`exp(-2π*im*e_quasi[i]*(t-t0)) * basis.modes(t)[i]`.
Periodicity of `H` is assumed. Solver options may be supplied directly or via
`propagator_kwargs`; direct keywords take precedence.
"""
function get_floquet_basis(H::AbstractQuantumObject, T::Real; t0::Real=0,
    propagator_kwargs=Dict{Symbol,Any}(), kwargs...)
    _floquet_period(T)
    isfinite(t0) || throw(ArgumentError("Reference time must be finite"))
    options = merge((; propagator_kwargs...), (; kwargs...))
    U = _floquet_propagator(H, T, t0; options...)
    values, states = eigenstates(U)
    energies = -angle.(values) ./ (2pi * T)
    order = sortperm(energies)
    energies, states = energies[order], states[order]
    modes = t -> propagate_floquet_modes(states, H, t, T;
        t0, e_quasi=energies, propagator_kwargs=options)
    return FloquetBasis(energies, modes, T, t0)
end
