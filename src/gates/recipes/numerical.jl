"""
    make_hamiltonian(prepared; scalar=1, remove_vacuum=false, herm_check=1e-9)

Compile the captured Hamiltonian, preserving `E/h` units with the default
`scalar=1`. Use `scalar=2pi` for angular-frequency units. The scalar and fixed
coefficient factors are folded into new matrices before `herm_check` cleanup.
Returns a QuantumObject without controls, otherwise QobjEvo. Controlled scalar
coefficients read selected canonical scalar keys from runtime `p`; unselected
pulse fields retain captured defaults. All
uncontrolled parameters keep their resolved build values. Fixed matrices and
scalar expression structure are compiled once.
`remove_vacuum=true` uses the captured model's constant reference vacuum energy;
added gate terms and pulse values do not redefine that reference. Coefficients retain
explicit adjoint terms; hermiticity cleanup follows `make_hamiltonian(model)`.
"""
function make_hamiltonian(prepared::PreparedGate; free_parameters=(), scalar::Number=1, remove_vacuum::Bool=false, herm_check=1e-9)
    return _compile_gate(prepared.model,prepared.recipe,prepared.terms,prepared.parameters;free_parameters,scalar,remove_vacuum,herm_check)
end

function _compile_gate(model,gate,extra,values;scalar::Number=1,remove_vacuum::Bool=false,herm_check=1e-9,runtime::Bool=false,free_parameters=())
    free = Set{Num}(c.target for c in gate.controls)
    bindings = runtime ? nothing : _compile_gate_bindings(gate,values,free_parameters)
    expressions = vcat([model.hamiltonian],extra)
    required = union((refs(expression,:param) for expression in expressions)...)
    issubset(free,required) || throw(ArgumentError("control target must appear in endpoint Hamiltonian or gate terms"))
    missing = setdiff(required,Set(keys(values)))
    isempty(missing) || throw(ArgumentError("missing reference values for parameters: $(join(string.(collect(missing)), ", "))"))
    term_list = reduce(vcat,physical_terms.(expressions);init=Tuple{Tuple,Sym}[])
    return compile_hamiltonian(term_list,model.operator_products,values,free,model_identity(model);
        scalar,herm_check,coefficient_bindings=bindings,
        energy_shift=remove_vacuum ? reference_vacuum_energy(model) : 0)
end

function make_hamiltonian(device::QuantumDevice,gate::GateRecipe;free_parameters=(),parameters=Dict{Num,Number}(),dressed_threads=nothing,tracking_threads=missing,kwargs...)
    prepared=build(device,gate;parameters,dressed_threads,tracking_threads)
    return make_hamiltonian(prepared;free_parameters,kwargs...)
end

function _runtime_gate_paths(gate,paths)
    paths isa AbstractString && throw(ArgumentError("free_parameters must be an ordered collection of scalar paths"))
    canonical=String[]
    for path in paths
        path isa AbstractString || throw(ArgumentError("runtime pulse paths must be strings"))
        index,name,keys=_gate_control_path(gate.controls,path)
        pulse=gate.controls[index].pulse
        _derived_pulse_path(pulse,keys) && throw(ArgumentError("derived duration cannot be selected: $name"))
        available=calibration_parameters(pulse)
        haskey(available,join(keys,'/')) || throw(ArgumentError("not a numeric scalar pulse field: $name"))
        _rebuild_only_parameter(pulse,keys) && throw(ArgumentError("$name requires rebuilding its envelope; use calibration candidate preparation or setpath"))
        name in canonical && throw(ArgumentError("duplicate runtime path or alias: $name"))
        push!(canonical,name)
    end
    return Tuple(canonical)
end
function _compile_gate_bindings(gate,values,paths)
    selected=_runtime_gate_paths(gate,paths)
    return Dict(c.target=>let
        recipe=_resolve_pulse(c.pulse.recipe,values)
        waveform=_compile_waveform(recipe,selected;base=string(c.name))
        (p,t) -> waveform(p,t,t)
    end for c in gate.controls)
end
function _prepare_gate_runtime(gate,reference,paths,params;canonical=nothing,default_trees=nothing)
    selected=isnothing(canonical) ? _runtime_gate_paths(gate,paths) : canonical
    runtime=Dict{String,Any}()
    for key in selected
        params!==nothing && haskey(params,key) || throw(ArgumentError("missing runtime pulse parameter $key"))
        value=params[key]
        value isa Number && !(value isa Bool) && concrete_parameter_value(value) || throw(ArgumentError("runtime pulse parameter $key must be a finite numeric scalar"))
        runtime[key]=value
    end
    runtime=_concrete_runtime_dict(runtime)
    params===nothing || all(k -> k in selected,keys(params)) || throw(ArgumentError("unknown runtime pulse parameter key"))
    resolved=map(Tuple(enumerate(gate.controls))) do (i,control)
        r=control.pulse.recipe
        defaults=isnothing(default_trees) ? _resolve_reference_tree(r,_parameter_tree(r),reference) : default_trees[i]
        _resolve_physical_parameters(r,defaults,runtime;base=string(control.name))
    end
    duration=gate.duration_policy===:auto ? maximum(v.delay+pulse_duration(c.pulse,v) for (c,v) in zip(gate.controls,resolved)) : gate.duration
    all(v.delay+pulse_duration(c.pulse,v)<=duration for (c,v) in zip(gate.controls,resolved)) || throw(ArgumentError("duration cuts off a runtime control"))
    stops=sort!(unique(vcat(Real[0,duration],[t for (c,v) in zip(gate.controls,resolved) for t in pulse_tstops(c.pulse,v)])))
    return (;runtime,resolved,duration,stops)
end
# Materialization is confined to returned provenance, after numerical evaluation.
function _effective_runtime_recipe(gate,resolved)
    controls=map(zip(gate.controls,resolved)) do (control,tree)
        entries=_numeric_parameter_dict(tree)
        changes=[(_path_keys(path),value) for (path,value) in entries if !_derived_pulse_path(control.pulse,split(path,'/')) && !_structural_pulse_parameter(control.pulse,split(path,'/'))]
        Control(control.name,control.target,Pulse(_pulse_update_tree(control.pulse.recipe,changes)))
    end
    return _physical_gate(gate.model_recipe,controls,gate.terms,gate.duration,gate.duration_policy,gate.metadata)
end

function _compile_gate_preparation(gate,reference,paths)
    canonical=_runtime_gate_paths(gate,paths)
    default_trees=map(gate.controls) do control
        r=control.pulse.recipe
        _resolve_reference_tree(r,_parameter_tree(r),reference)
    end
    return params -> _prepare_gate_runtime(gate,reference,canonical,params;canonical,default_trees)
end
