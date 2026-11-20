function _floquet_period(period)
    period isa Real && !(period isa Bool) && isfinite(period) && period>0 || throw(ArgumentError("Floquet period must be finite and positive"))
    T=Float64(period)
    isfinite(T) && T>0 && isfinite(inv(T)) || throw(ArgumentError("period must be representable with a finite inverse frequency"))
    return T
end
function _floquet_inputs(H,period,t0,params)
    T=_floquet_period(period)
    t0 isa Real && isfinite(t0) || throw(ArgumentError("reference time must be finite"))
    H isa Union{QuantumObject,QobjEvo} || throw(ArgumentError("Floquet Hamiltonian must be a QuantumObject or QobjEvo"))
    reference=H isa QobjEvo ? H(params,t0) : H
    isoper(reference) && size(reference,1)==size(reference,2) || throw(ArgumentError("Floquet Hamiltonian must be a square operator"))
    all(isfinite,reference.data) || throw(ArgumentError("Hamiltonian contains nonfinite values"))
    return reference,T,Float64(t0)
end
function _floquet_solver_options(settings;final_time)
    settings isa NamedTuple || throw(ArgumentError("evolution_settings must be a NamedTuple"))
    reserved=(:e_ops,:save_idxs,:saveat,:save_start,:save_everystep,:callback,:params)
    any(k->haskey(settings,k),reserved) && throw(ArgumentError("evolution settings cannot override full-propagator output or params"))
    get(settings,:save_end,true)===true || throw(ArgumentError("Floquet propagation requires save_end=true"))
    return merge((;progress_bar=false),settings,(;saveat=[final_time],save_start=false,save_end=true,save_everystep=false))
end
function _floquet_tolerance(value,name)
    value isa Real && isfinite(value) && value>=0 || throw(ArgumentError("$name must be finite and nonnegative"))
    return Float64(value)
end
