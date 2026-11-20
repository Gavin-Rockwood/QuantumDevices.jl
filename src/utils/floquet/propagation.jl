"""
    propagate_floquet_modes(basis, H, t; params=nothing, evolution_settings=(;))

Return periodic modes at absolute time t. Supply the same periodic drive and
parameters explicitly; a saved basis contains neither. No physical dynamical
phase over whole periods is included; use propagate_floquet_states for that.
"""
function propagate_floquet_modes(basis::FloquetBasis,H,t;params=nothing,evolution_settings=(;))
    t isa Real && isfinite(t) || throw(ArgumentError("Floquet sampling time must be finite"))
    reference,_,_=_floquet_inputs(H,basis.period,basis.t0,params)
    reference.dimensions.to==first(basis.states).dimensions.to || throw(DimensionMismatch("drive and Floquet reference dimensions differ"))
    elapsed=mod(t-basis.t0,basis.period)
    iszero(elapsed) && return deepcopy(basis.states)
    settings=merge(basis.settings.evolution,evolution_settings)
    U,_=_floquet_propagator(H,elapsed,basis.t0;params,evolution_settings=settings)
    modes=[cis(2pi*basis.e_quasi[i]*elapsed)*(U*state) for (i,state) in enumerate(basis.states)]
    vectors=hcat((state.data for state in modes)...)
    norm(vectors'*vectors-I)/sqrt(length(modes))<=basis.settings.checks.mode_atol || error("Micromotion Gram defect exceeds tolerance; tighten evolution tolerances")
    return modes
end
"""Return physical Floquet solutions, including exp(-2πim*e*(t-t0)), at absolute time t."""
function propagate_floquet_states(basis::FloquetBasis,H,t;kwargs...)
    modes=propagate_floquet_modes(basis,H,t;kwargs...)
    return [cis(-2pi*basis.e_quasi[i]*(t-basis.t0))*state for (i,state) in enumerate(modes)]
end
