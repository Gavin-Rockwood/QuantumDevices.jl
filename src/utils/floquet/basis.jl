_fold_quasienergy(energy,period)=mod(energy+0.5/period,1/period)-0.5/period
_circular_gap(a,b,frequency)=abs(mod(a-b+frequency/2,frequency)-frequency/2)
function _circular_mode_diagnostics(energies,period;atol,rtol)
    frequency=1/period
    nearest=fill(Inf,length(energies));second=copy(nearest)
    for i in eachindex(energies)
        distances=sort!([_circular_gap(energies[i],energies[j],frequency) for j in eachindex(energies) if j!=i])
        isempty(distances) || (nearest[i]=distances[1])
        length(distances)>1 && (second[i]=distances[2])
    end
    return (;nearest_gap=nearest,second_gap=second,floquet_degenerate=nearest .<= atol+rtol*frequency)
end
# Propagation owns its identity/scratch and applies the cycles-to-angular factor once.
function _floquet_propagator(H,duration,t0;params=nothing,evolution_settings=(;))
    reference,_,_=_floquet_inputs(H,max(duration,eps(Float64)),t0,params)
    initial=QuantumObject(Matrix{ComplexF64}(I,size(reference,1),size(reference,1));dims=reference.dimensions)
    iszero(duration) && return initial,(;retcode=SciMLBase.ReturnCode.Success)
    isfinite(t0+duration) && t0+duration>t0 || throw(ArgumentError("period endpoint must be finite and distinguishable from t0"))
    options=_floquet_solver_options(evolution_settings;final_time=t0+duration)
    solution=sesolve(2pi*H,initial,[t0,t0+duration];params,options...)
    SciMLBase.successful_retcode(solution.retcode) || error("Floquet propagation failed: $(solution.retcode)")
    !isempty(solution.states) && last(solution.times_states)==t0+duration || error("Floquet propagation did not save the complete final operator")
    return last(solution.states),(;retcode=solution.retcode,alg=solution.alg,abstol=solution.abstol,reltol=solution.reltol)
end
function _floquet_eigensystem(U,period;unitarity_atol,mode_atol)
    matrix=Matrix(U.data);n=size(matrix,1);scale=sqrt(n)
    unitarity=norm(matrix'*matrix-I)/scale
    decomposition=schur(matrix)
    values=diag(decomposition.T);vectors=decomposition.Z
    radius=maximum(abs.(abs.(values).-1))
    gram=norm(vectors'*vectors-I)/scale
    residual=norm(matrix*vectors-vectors*Diagonal(values))/scale
    offdiagonal=norm(decomposition.T-Diagonal(values))/scale
    all(isfinite,(unitarity,radius,gram,residual,offdiagonal)) || error("Nonfinite Floquet decomposition")
    max(unitarity,radius)<=unitarity_atol || error("Floquet propagator is not unitary within tolerance; tighten evolution tolerances")
    max(gram,residual,offdiagonal)<=mode_atol || error("Floquet mode residual/orthogonality exceeds tolerance; no matrix repair was applied")
    energies=_fold_quasienergy.(-angle.(values)./(2pi*period),period)
    order=sortperm(energies;alg=Base.Sort.MergeSort)
    dimensions=QuantumToolbox.Dimensions(U.dimensions.to,QuantumToolbox.Space(1))
    states=[QuantumObject(copy(vectors[:,i]);type=Ket(),dims=dimensions) for i in order]
    return energies[order],states,(;unitarity,radius,gram,residual,schur_offdiagonal=offdiagonal)
end
"""
    get_floquet_basis(H, period; t0=0, params=nothing, evolution_settings=(;), ...)

Compute reference modes and quasienergies in [-1/(2period),1/(2period)). H uses
cycles/time and is assumed periodic. Results contain data, never a drive closure.
Nonunitary propagation is diagnosed and rejected, not silently repaired.
"""
function get_floquet_basis(H,period;t0=0,params=nothing,evolution_settings=(;),unitarity_atol=1e-7,mode_atol=1e-7,degeneracy_atol=1e-9,degeneracy_rtol=1e-10)
    _,T,origin=_floquet_inputs(H,period,t0,params)
    checks=(;unitarity_atol=_floquet_tolerance(unitarity_atol,"unitarity_atol"),mode_atol=_floquet_tolerance(mode_atol,"mode_atol"),
        degeneracy_atol=_floquet_tolerance(degeneracy_atol,"degeneracy_atol"),degeneracy_rtol=_floquet_tolerance(degeneracy_rtol,"degeneracy_rtol"))
    U,solver=_floquet_propagator(H,T,origin;params,evolution_settings)
    energies,states,defects=_floquet_eigensystem(U,T;unitarity_atol=checks.unitarity_atol,mode_atol=checks.mode_atol)
    circular=_circular_mode_diagnostics(energies,T;atol=checks.degeneracy_atol,rtol=checks.degeneracy_rtol)
    return FloquetBasis(energies,states,T,origin,merge(defects,circular),(;evolution=deepcopy(evolution_settings),solver,checks))
end
