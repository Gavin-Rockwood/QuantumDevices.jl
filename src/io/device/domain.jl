import DifferentialEquations, ADTypes
const _IO_DOMAIN_TYPES = Dict{String,Any}(
    "OperatorGenerator"=>OperatorGenerator,"Input"=>Input,"InputBinding"=>InputBinding,
    "ResolvedRecipe"=>ResolvedRecipe,"Control"=>Control,"LabelDiagnostic"=>LabelDiagnostic,
    "DressedStates"=>DressedStates,"PreparedGate"=>PreparedGate,"GateEvolution"=>GateEvolution)
function _fields_record(x,w,path;omit=())
    return Dict(String(k)=>_encode_record(getfield(x,k),w,"$path/$k") for k in fieldnames(typeof(x)) if !(k in omit))
end
function _domain_record(x,w,path)
    x isa _AnalysisArtifact && return _analysis_record(x,w,path)
    x isa ADTypes.AutoFiniteDiff && return Dict("tag"=>"auto_finite_diff","fields"=>_fields_record(x,w,path))
    x isa Val && return Dict("tag"=>"val","value"=>_encode_record(first(typeof(x).parameters),w,path))
    for id in (:Vern6,:Vern7,:Vern8,:Vern9,:Tsit5,:DP5)
        isdefined(DifferentialEquations,id) || continue
        T=getproperty(DifferentialEquations,id)
        if x isa T
            default=T()
            for field in fieldnames(typeof(x))
                field===:lazy && continue
                isequal(getfield(x,field),getfield(default,field)) || return _custom_record(x,w,path)
            end
            lazy=hasfield(typeof(x),:lazy) ? (x.lazy isa Bool ? x.lazy : first(typeof(x.lazy).parameters)) : nothing
            return Dict("tag"=>"ode_algorithm","id"=>String(id),"lazy"=>lazy)
        end
    end
    if x isa Component
        f=Dict("name"=>_encode_record(x.name,w,path),"variables"=>_encode_record(x.variables,w,path),"parameters"=>_encode_record(x.parameters,w,path),
          "operators"=>_encode_record(Dict(x.operators),w,path),"operator_generators"=>_encode_record(Dict(x.operator_generators),w,path),
          "hamiltonian"=>_encode_record(x.hamiltonian,w,path),"type"=>_encode_record(x.type,w,path),"label"=>_encode_record(x.label,w,path),"metadata"=>_encode_record(x.metadata,w,path))
        return Dict("tag"=>"domain","kind"=>"Component","fields"=>f)
    end
    if x isa GateRecipe
        return Dict("tag"=>"domain","kind"=>"GateRecipe","fields"=>_fields_record(x,w,path))
    end
    if x isa CalibrationResult
        summary=(;u=copy(collect(x.solution.u)),retcode=x.solution.retcode)
        hasproperty(x.solution,:objective) && (summary=merge(summary,(;objective=x.solution.objective)))
        f=_fields_record(x,w,path;omit=(:solution,))
        f["solution"]=_encode_record(summary,w,path*"/solution")
        id="result$(length(w.native)+1)"
        w.native[id]=x.solution
        return Dict("tag"=>"domain","kind"=>"CalibrationResult","fields"=>f,"native_id"=>id)
    end
    if x isa TrackingResult
        return Dict("tag"=>"domain","kind"=>"TrackingResult","fields"=>_fields_record(x,w,path;omit=(:_lookup,)))
    end
    if x isa TrackedData
        owner=x.tracking;storage=getfield(x,:_storage)
        selector=nothing
        for field in (:_states,:_indices,:_diagnostics)
            storage===getfield(owner,field) && (selector=(field,))
        end
        for field in (:_other,:_accuracy)
            group=getfield(owner,field)
            isnothing(group) && continue
            for (key,value) in pairs(group)
                storage===value && (selector=(field,key))
            end
        end
        return Dict("tag"=>"tracked_data","tracking"=>_encode_record(owner,w,path*"/tracking"),
            "selector"=>_encode_record(selector,w,path),"storage"=>isnothing(selector) ? _encode_record(storage,w,path) : nothing)
    end
    if x isa Union{Interaction,NumericalModel,QuantumGate}
        return Dict("tag"=>"domain","kind"=>string(nameof(typeof(x))),"fields"=>_fields_record(x,w,path))
    end
    name=string(nameof(typeof(x)))
    if haskey(_IO_DOMAIN_TYPES,name) && x isa _IO_DOMAIN_TYPES[name]
        fields=_fields_record(x,w,path)
        x isa InputBinding && (fields["scope"]=Dict("tag"=>"scope","id"=>_scope(w,x.scope)))
        return Dict("tag"=>"domain","kind"=>name,"fields"=>fields)
    end
    if x isa _PulseObject
        id=findfirst(T->x isa T,_PULSE_IMPLEMENTATIONS)
        isnothing(id) && return _custom_record(x,w,path)
        fields=Dict(String(k)=>_encode_record(getproperty(x,k),w,"$path/$k") for k in _pulse_storage_fields(x))
        return Dict("tag"=>"pulse","version"=>2,"id"=>String(id),"fields"=>fields)
    end
    return _custom_record(x,w,path)
end
function _restore_domain(d,r,path)
    tag=d["tag"]
    tag=="analysis" && return _restore_analysis(d,r,path)
    if tag=="ode_algorithm"
        id=Symbol(d["id"]);id in (:Vern6,:Vern7,:Vern8,:Vern9,:Tsit5,:DP5) || _ioerror(path,"unknown algorithm")
        T=getproperty(DifferentialEquations,id)
        return isnothing(d["lazy"]) ? T() : T(;lazy=Bool(d["lazy"]))
    end
    tag=="val" && return Val(_decode_record(d["value"],r,path))
    if tag=="auto_finite_diff"
        f=(; (Symbol(k)=>_decode_record(v,r,path) for (k,v) in d["fields"])...)
        return ADTypes.AutoFiniteDiff(;f...)
    end
    if tag=="tracked_data"
        owner=_decode_record(d["tracking"],r,path);selector=_decode_record(d["selector"],r,path)
        storage=isnothing(selector) ? _decode_record(d["storage"],r,path) : length(selector)==1 ? getfield(owner,selector[1]) : getproperty(getfield(owner,selector[1]),selector[2])
        return TrackedData(storage,owner)
    end
    tag=="scope" && return _unScope(r,d["id"])
    if tag=="pulse"
        f=(; (Symbol(k)=>_decode_record(v,r,"$path/$k") for (k,v) in d["fields"])...)
        try
            return _restore_pulse_fields(Symbol(d["id"]),get(d,"version",1),f)
        catch error
            error isa InterruptException && rethrow()
            _ioerror(path,sprint(showerror,error))
        end
    end
    tag=="domain" || _ioerror(path,"unknown record tag $tag")
    kind=d["kind"]
    types=merge(_IO_DOMAIN_TYPES,Dict("Component"=>Component,"Interaction"=>Interaction,"GateRecipe"=>GateRecipe,
        "QuantumGate"=>QuantumGate,"NumericalModel"=>NumericalModel,"TrackingResult"=>TrackingResult,"TrackedData"=>TrackedData,"CalibrationResult"=>CalibrationResult))
    haskey(types,kind) || _ioerror(path,"unknown domain type $kind")
    expected=Set(String.(fieldnames(types[kind])))
    kind=="TrackingResult" && delete!(expected,"_lookup")
    Set(keys(d["fields"]))==expected || _ioerror(path,"domain field schema mismatch for $kind")
    f=(; (Symbol(k)=>_decode_record(v,r,"$path/$k") for (k,v) in d["fields"])...)
    kind=="Component" && return Component(f.name;variables=f.variables,parameters=f.parameters,operators=f.operators,operator_generators=f.operator_generators,hamiltonian=f.hamiltonian,type=f.type,label=f.label,metadata=f.metadata)
    kind=="Interaction" && return Interaction(f.name,f.expression,f.variables,f.parameters,f.metadata)
    kind=="GateRecipe" && return _physical_gate(f.model_recipe,f.controls,f.terms,f.duration,f.duration_policy,f.metadata)
    kind=="QuantumGate" && return QuantumGate(f.name,_io_working(f.recipe),f.calibration,f.metadata,Dict{Symbol,Any}(f.artifacts))
    kind=="PreparedGate" && return PreparedGate(f.recipe,f.model,Sym[f.terms...],Dict{Num,Number}(f.parameters),f.context)
    kind=="NumericalModel" && return NumericalModel(f.recipe,f.H,f.eigensystem,Dict{OperatorRef,QuantumObject}(f.operators),f.hamiltonian,Set{Symbol}(f.owners),Dict{Tuple,QuantumObject}(f.operator_products))
    kind=="TrackingResult" && return TrackingResult(f.keys,f.saved_steps,f.total_steps,f._states,f._other,f._indices,f._accuracy,f._diagnostics;tracked=f.tracked)
    kind=="TrackedData" && return TrackedData(f._storage,f.tracking)
    kind=="CalibrationResult" && return CalibrationResult(f.recipe,get(r.native,get(d,"native_id",""),f.solution),f.coordinates,Tuple(String.(f.paths)),f.parameterization,f.fitted_values,f.state_keys,f.states,f.target,f.virtual_z,f.settings,f.provenance,f.validation)
    kind=="InputBinding" && return InputBinding(f.name,f.recipe,_io_owned_vector(f.owners),f.scope,_io_owned_map(f.operator_names))
    kind=="ResolvedRecipe" && return ResolvedRecipe(f.definition,_io_owned_map(f.parameters),_io_owned_vector(f.inputs),Val(:stored))
    kind=="Input" && return Input(f.name,f.recipe,_io_owned_vector(f.states))
    haskey(_IO_DOMAIN_TYPES,kind) || _ioerror(path,"unknown domain type $kind")
    T=_IO_DOMAIN_TYPES[kind]
    return T((getproperty(f,k) for k in fieldnames(T))...)
end
