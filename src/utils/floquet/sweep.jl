function _floquet_sample(H,period,time;t0,params,evolution_settings,floquet_kwargs)
    basis=get_floquet_basis(H,period;t0,params,evolution_settings,floquet_kwargs...)
    modes=propagate_floquet_modes(basis,H,time;params)
    validate_kets(modes) # Tracking owns its existing normalization contract; do not renormalize.
    d=basis.diagnostics
    defects=(;unitarity=d.unitarity,radius=d.radius,gram=d.gram,residual=d.residual,schur_offdiagonal=d.schur_offdiagonal)
    return (;modes,energies=basis.e_quasi,nearest=d.nearest_gap,second=d.second_gap,degenerate=d.floquet_degenerate,defects)
end
"""
    floquet_sweep(H_func, points, periods; seed=nothing, keys=nothing, save_steps=:all, ...)

Compute periodic modes in a common construction basis and match them with
track_states. Returned selected histories have key rows and saved-step columns.
Factories run once per point, without caching. Threaded factories must return
independent mutable workspaces. All candidate histories are transiently retained.
"""
function floquet_sweep(H_func,points,periods;t0=0,sample_times=nothing,params=nothing,seed=nothing,keys=nothing,save_steps=:all,
    floquet_threads=1,tracking_kwargs=(;),evolution_settings=(;),floquet_kwargs=(;))
    samples=collect(points);n=length(samples)
    n>0 || throw(ArgumentError("Floquet sweep requires nonempty sampling points"))
    T=periods isa Real ? fill(_floquet_period(periods),n) : _floquet_period.(collect(periods))
    length(T)==n || throw(DimensionMismatch("supply one period per point"))
    times=isnothing(sample_times) ? fill(Float64(t0),n) : Float64.(collect(sample_times))
    length(times)==n && all(isfinite,times) || throw(ArgumentError("supply one finite sampling time per point"))
    floquet_threads isa Int && floquet_threads>0 || throw(ArgumentError("floquet_threads must be a positive Int"))
    tracking_kwargs isa NamedTuple && floquet_kwargs isa NamedTuple || throw(ArgumentError("tracking_kwargs and floquet_kwargs must be NamedTuples"))
    any(k->haskey(tracking_kwargs,k),(:seed,:keys,:other,:save_steps,:transports,:seed_transport)) && throw(ArgumentError("tracking_kwargs cannot override sweep references or saved histories"))
    any(k->haskey(floquet_kwargs,k),(:t0,:params,:evolution_settings)) && throw(ArgumentError("floquet_kwargs cannot override time/params/evolution settings"))
    results=Vector{Any}(undef,n)
    function sample!(i)
        try
            results[i]=_floquet_sample(H_func(samples[i]),T[i],times[i];t0,params,evolution_settings,floquet_kwargs)
        catch error
            throw(ErrorException("Floquet sample $i failed: $(sprint(showerror,error))"))
        end
    end
    workers=min(floquet_threads,Threads.nthreads(:default),n)
    if workers==1
        foreach(sample!,1:n)
    else
        @sync for worker in 1:workers
            Threads.@spawn for i in worker:workers:n;sample!(i);end
        end
    end
    history=[r.modes for r in results]
    other=(;quasienergies=[r.energies for r in results],nearest_gap=[r.nearest for r in results],second_gap=[r.second for r in results],floquet_degenerate=[r.degenerate for r in results])
    tracking=track_states(history;seed,keys,save_steps,other,tracking_kwargs...)
    settings=(;t0=Float64(t0),evolution=deepcopy(evolution_settings),floquet=deepcopy(floquet_kwargs),tracking=deepcopy(tracking_kwargs),floquet_threads)
    return FloquetSweep(samples,T,times,tracking,[r.defects for r in results],settings)
end
