# Preserve tensor and Liouville structure, rather than inferring it from matrix size.
function _space_record(space,w,path)
    name=string(nameof(typeof(space)))
    name=="Space" && return Dict("kind"=>name,"size"=>space.size)
    name=="TensorSpace" && return Dict("kind"=>name,"spaces"=>[_space_record(s,w,path) for s in space.spaces])
    name=="LiouvilleSpace" && return Dict("kind"=>name,"to"=>_space_record(space.op_dims.to,w,path),"from"=>_space_record(space.op_dims.from,w,path))
    return Dict("kind"=>"custom","record"=>_custom_record(space,w,path))
end
function _restore_space(d,r,path)
    d["kind"]=="Space" && return QuantumToolbox.Space(Int(d["size"]))
    d["kind"]=="TensorSpace" && return QuantumToolbox.TensorSpace(Tuple(_restore_space(s,r,path) for s in d["spaces"]))
    d["kind"]=="LiouvilleSpace" && return QuantumToolbox.LiouvilleSpace(QuantumToolbox.Dimensions(_restore_space(d["to"],r,path),_restore_space(d["from"],r,path)))
    d["kind"]=="custom" && return _decode_record(d["record"],r,path)
    _ioerror(path,"unknown quantum space")
end
function _quantum_record(q,w,path)
    return Dict("tag"=>"quantum","role"=>string(nameof(typeof(q.type))),
        "to"=>_space_record(q.dimensions.to,w,path*"/to"),"from"=>_space_record(q.dimensions.from,w,path*"/from"),"data"=>_encode_record(q.data,w,path*"/data"))
end
function _restore_quantum(d,r,path)
    to=_restore_space(d["to"],r,path);from=_restore_space(d["from"],r,path)
    data=_decode_record(d["data"],r,path*"/data")
    role=get(Dict("Ket"=>Ket,"Bra"=>Bra,"Operator"=>Operator,"SuperOperator"=>SuperOperator,"OperatorKet"=>OperatorKet,"OperatorBra"=>OperatorBra),d["role"],nothing)
    isnothing(role) && _ioerror(path,"unsupported quantum role")
    return QuantumObject(data;type=role(),dims=QuantumToolbox.Dimensions(to,from))
end
