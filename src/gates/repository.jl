# Gate paths belong to the repository, independently of a gate's own name.
function _gate_path(path)
    result = path isa Symbol ? (path,) : path isa Tuple ? path : ()
    !isempty(result) && all(x -> x isa Symbol,result) || throw(ArgumentError("gate path must be a Symbol or nonempty tuple of Symbols"))
    return result
end
function _validate_gate_repository(repository, active=Base.IdSet())
    repository in active && throw(ArgumentError("cyclic gate repository"))
    push!(active,repository)
    for (key,value) in repository
        key isa Symbol || throw(ArgumentError("gate repository keys must be Symbols"))
        value isa QuantumGate || value isa Dict || throw(ArgumentError("gate repository values must be gates or dictionaries"))
        value isa Dict && _validate_gate_repository(value,active)
    end
    delete!(active,repository)
    return repository
end
"""Return the gate at a Symbol or nested tuple path; groups and missing paths error."""
function get_gate(device::QuantumDevice,path)
    node = device.gates
    for key in _gate_path(path)
        node isa AbstractDict && haskey(node,key) || throw(ArgumentError("unknown gate path: $path"))
        node = node[key]
    end
    node isa QuantumGate || throw(ArgumentError("gate path names a group: $path"))
    return node
end
function _gate_parent!(repository,path)
    node=repository
    for key in path[1:end-1]
        haskey(node,key) || (node[key]=Dict{Symbol,Union{Dict,QuantumGate}}())
        node[key] isa Dict || throw(ArgumentError("gate parent is a leaf: $key"))
        node=node[key]
    end
    return node
end
# Copy dictionary structure only, retaining leaf identity during preflight.
_gate_tree_copy(repository) = Dict{Symbol,Union{Dict,QuantumGate}}(k=>(v isa Dict ? _gate_tree_copy(v) : v) for (k,v) in repository)
function _insert_gate!(repository,path,gate;replace=false)
    parent=_gate_parent!(repository,path)
    existing=haskey(parent,last(path))
    replace ? (existing && parent[last(path)] isa QuantumGate || throw(ArgumentError("replacement requires an existing gate"))) :
        (!existing || throw(ArgumentError("gate destination already exists")))
    parent[last(path)]=gate
    return repository
end
"""Register a gate at a unique path, creating groups; `add_gate!(device,gate)` uses its name."""
function add_gate!(device::QuantumDevice,path,gate::QuantumGate)
    _validate_gate_repository(device.gates)
    keys=_gate_path(path)
    _insert_gate!(_gate_tree_copy(device.gates),keys,gate)
    _insert_gate!(device.gates,keys,gate)
    return gate
end
add_gate!(device::QuantumDevice,gate::QuantumGate)=add_gate!(device,gate.name,gate)
"""Explicitly replace an existing gate leaf without renaming the supplied gate."""
function replace_gate!(device::QuantumDevice,path,gate::QuantumGate)
    _validate_gate_repository(device.gates)
    keys=_gate_path(path)
    _insert_gate!(_gate_tree_copy(device.gates),keys,gate;replace=true)
    _insert_gate!(device.gates,keys,gate;replace=true)
    return gate
end
"""Atomically move a gate leaf between repository paths, preserving object identity and snapshots."""
function move_gate!(device::QuantumDevice,source,destination)
    _validate_gate_repository(device.gates)
    src,dst=_gate_path(source),_gate_path(destination)
    gate=get_gate(device,src)
    src==dst && return gate
    prospective=_gate_tree_copy(device.gates)
    delete!(_gate_parent!(prospective,src),last(src))
    _insert_gate!(prospective,dst,gate)
    # All dictionary/type/path validation has completed before touching the original.
    old=_gate_tree_copy(device.gates)
    try
        empty!(device.gates);merge!(device.gates,prospective)
    catch
        empty!(device.gates);merge!(device.gates,old);rethrow()
    end
    return gate
end
"""Attach a named saved evolution, tracking result, model or calibration; never runs a simulation."""
function attach_artifact!(gate::QuantumGate,name::Symbol,artifact;replace=false)
    artifact isa Union{GateEvolution,TrackingResult,DressedStates,NumericalModel,CalibrationResult,_AnalysisArtifact} ||
        throw(ArgumentError("unsupported gate artifact $(typeof(artifact))"))
    !haskey(gate.artifacts,name) || replace || throw(ArgumentError("artifact already exists: $name"))
    gate.artifacts[name]=artifact
    return artifact
end
