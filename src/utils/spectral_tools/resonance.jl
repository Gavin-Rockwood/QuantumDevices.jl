function _resonance_references(references,keys)
    references isa Union{AbstractVector,Tuple} || throw(ArgumentError("reference states must be an ordered ket collection"))
    states=collect(references)
    length(states)==2 || throw(ArgumentError("select exactly two reference states"))
    validate_kets(states)
    abs(dot(states[1],states[2]))<=1e-8 || throw(ArgumentError("reference states must be orthogonal"))
    names=isnothing(keys) ? ((0,),(1,)) : Tuple(keys)
    length(names)==2 && allunique(names) || throw(ArgumentError("supply two unique reference keys"))
    return deepcopy(states),names
end
function _resonance_result(sweep;fit_kwargs,reference_frequency,photon_order,metadata)
    energies=sweep.tracking.other.quasienergies
    gaps=[_circular_gap(energies[1,j],energies[2,j],sweep.points[j]) for j in eachindex(sweep.points)]
    fit=fit_avoided_crossing(sweep.points,gaps;fit_kwargs)
    fit.minimum_gap>sqrt(eps(Float64))*maximum(gaps) || error("Fitted resonance gap is unresolved")
    # A nearby third mode is a warning sign, not proof of isolated two-state dynamics.
    center=argmin(abs.(sweep.points .-fit.center))
    neighborhood=max(1,center-1):min(length(gaps),center+1)
    competing=any(j->minimum(sweep.tracking.other.second_gap[:,j])<=3gaps[j],neighborhood)
    diagnostics=(;potential_competing_modes=competing,low_overlap=any(sweep.tracking.accuracy.low_overlap),
        circular_degeneracy=any(sweep.tracking.other.floquet_degenerate),competing_gap_factor=3.,competing_samples=Tuple(neighborhood),two_state_assumption=:unchecked)
    return ResonanceResult(sweep,fit,reference_frequency,photon_order,metadata,diagnostics)
end
"""
    find_resonance(H_func, frequencies, reference_states; keys=nothing, ...)

Fit the circular gap of two tracked modes for drives with fundamental period 1/f.
Requires a bracketed resolved crossing. Residuals and overlaps diagnose quality;
they do not certify isolated two-state dynamics or a calibrated pulse duration.
"""
function find_resonance(H_func,frequencies,reference_states;keys=nothing,t0=0,params=nothing,photon_order=1,reference_frequency=nothing,
    floquet_threads=1,tracking_kwargs=(;),evolution_settings=(;),floquet_kwargs=(;),fit_kwargs=(;))
    freqs,_=_spectral_parameters(frequencies)
    all(>(0),freqs) || throw(ArgumentError("drive frequencies must be positive"))
    photon_order isa Integer && !(photon_order isa Bool) && photon_order>0 || throw(ArgumentError("photon_order must be a positive integer"))
    states,names=_resonance_references(reference_states,keys)
    sweep=floquet_sweep(H_func,freqs,1 ./freqs;t0,params,seed=states,keys=names,save_steps=:all,
        floquet_threads,tracking_kwargs,evolution_settings,floquet_kwargs)
    return _resonance_result(sweep;fit_kwargs,reference_frequency,photon_order,metadata=(;drive=:factory))
end
"""Scan H0 + amplitude*cos(2π*f*t+phase)*drive; phase uses absolute time, not t0."""
function find_resonance(H0::QuantumObject,drive::QuantumObject,amplitude,frequencies,reference_states;phase=0,kwargs...)
    isoper(H0) && isoper(drive) && H0.dimensions==drive.dimensions || throw(DimensionMismatch("H0 and drive must be compatible operators"))
    all(x->x isa Real && isfinite(x),(amplitude,phase)) || throw(ArgumentError("amplitude and phase must be finite real values"))
    norm(H0.data-H0.data')<=1e-9 && norm(drive.data-drive.data')<=1e-9 || throw(ArgumentError("H0 and cosine drive must be Hermitian"))
    factory=f->QobjEvo((copy(H0),(copy(drive),(p,t)->amplitude*cos(2pi*f*t+phase))))
    result=find_resonance(factory,frequencies,reference_states;kwargs...)
    return ResonanceResult(result.sweep,result.fit,result.reference_frequency,result.photon_order,(;drive=:cosine,amplitude,phase),result.diagnostics)
end
