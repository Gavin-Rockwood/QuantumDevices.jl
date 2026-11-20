"""Data-only reference Floquet modes and principal-zone quasienergies in cycles/time."""
struct FloquetBasis{E,S,D,K}
    e_quasi::E
    states::S
    period::Float64
    t0::Float64
    diagnostics::D
    settings::K
    function FloquetBasis(energies,states,period,t0,diagnostics,settings)
        _floquet_period(period)
        isfinite(t0) || throw(ArgumentError("reference time must be finite"))
        length(energies)==length(states)>0 || throw(DimensionMismatch("mode and energy counts differ"))
        all(e->isfinite(e) && -0.5/period<=e<0.5/period,energies) || throw(ArgumentError("quasienergies must lie in the half-open principal zone"))
        diagnostics isa NamedTuple && all(k->haskey(diagnostics,k),(:unitarity,:radius,:gram,:residual,:schur_offdiagonal,:nearest_gap,:second_gap,:floquet_degenerate)) || throw(ArgumentError("incomplete Floquet diagnostics"))
        all(k->length(getproperty(diagnostics,k))==length(states),(:nearest_gap,:second_gap,:floquet_degenerate)) || throw(DimensionMismatch("mode diagnostic counts differ"))
        all(k->getproperty(diagnostics,k) isa Real && isfinite(getproperty(diagnostics,k)) && getproperty(diagnostics,k)>=0,(:unitarity,:radius,:gram,:residual,:schur_offdiagonal)) || throw(ArgumentError("invalid Floquet defects"))
        all(x->x isa Real && !isnan(x) && x>=0,diagnostics.nearest_gap) && all(x->x isa Real && !isnan(x) && x>=0,diagnostics.second_gap) && all(x->x isa Bool,diagnostics.floquet_degenerate) || throw(ArgumentError("invalid circular diagnostics"))
        all(k->haskey(settings,k),(:evolution,:solver,:checks)) || throw(ArgumentError("incomplete Floquet settings"))
        all(k->haskey(settings.checks,k),(:mode_atol,:unitarity_atol,:degeneracy_atol,:degeneracy_rtol)) || throw(ArgumentError("incomplete check tolerances"))
        validate_kets(states)
        length(states)==length(first(states).data) || throw(DimensionMismatch("Floquet basis must span the reference space"))
        new{typeof(energies),typeof(states),typeof(diagnostics),typeof(settings)}(energies,states,Float64(period),Float64(t0),diagnostics,settings)
    end
end
"""Ordered periodic-drive samples, with selected trajectories owned by current TrackingResult."""
struct FloquetSweep{P,T,D,K}
    points::P
    periods::Vector{Float64}
    sample_times::Vector{Float64}
    tracking::T
    diagnostics::D
    settings::K
    function FloquetSweep(points,periods,times,tracking,diagnostics,settings)
        n=length(points)
        n>0 && length(periods)==length(times)==length(diagnostics)==tracking.total_steps==n || throw(DimensionMismatch("Floquet sweep sample counts differ"))
        foreach(_floquet_period,periods)
        all(isfinite,times) || throw(ArgumentError("sampling times must be finite"))
        all(k->haskey(getfield(tracking,:_other),k),(:quasienergies,:nearest_gap,:second_gap,:floquet_degenerate)) || throw(ArgumentError("Floquet tracking requires quasienergies"))
        stored=getfield(tracking,:_states);energies=getfield(tracking,:_other).quasienergies
        for (column,step) in enumerate(tracking.saved_steps)
            kets=ndims(stored)==1 ? stored : stored[:,column]
            validate_kets(kets)
            values=ndims(energies)==1 ? energies : energies[:,column]
            all(e->e isa Real && isfinite(e) && -0.5/periods[step]<=e<0.5/periods[step],values) || throw(ArgumentError("invalid saved quasienergies"))
        end
        new{typeof(points),typeof(tracking),typeof(diagnostics),typeof(settings)}(points,Float64.(periods),Float64.(times),tracking,diagnostics,settings)
    end
end
