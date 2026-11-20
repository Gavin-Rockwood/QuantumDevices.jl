import SymbolicUtils
const _IO_VARIABLE_METADATA = Dict(string(name)=>getproperty(Symbolics,name) for name in (:VariableSource,:VariableDefaultValue) if isdefined(Symbolics,name))
function _scalar_metadata(raw,w,path)
    metadata=SymbolicUtils.metadata(raw)
    isnothing(metadata) && return Dict{String,Any}()
    result=Dict{String,Any}()
    for (key,value) in metadata
        id=findfirst(isequal(key),_IO_VARIABLE_METADATA)
        if isnothing(id)
            key_record=_custom_record(key,w,path*"/metadata/key")
            result["custom:"*recipe_fingerprint(key_record)]=Dict("key"=>key_record,"value"=>_encode_record(value,w,path*"/metadata/value"))
        else
            result[id]=_encode_record(value,w,path*"/metadata/$id")
        end
    end
    return result
end
function _encode_scalar(x,w,path)
    raw=SymbolicUtils.unwrap_const(Symbolics.unwrap(x))
    if Symbolics.iscall(raw)
        operation=Symbolics.operation(raw)
        id=findfirst(isequal(operation),_IO_OPERATIONS)
        operation_record=isnothing(id) ? _custom_record(operation,w,path*"/operation") : Dict("tag"=>"operation","id"=>id)
        return Dict("tag"=>"scalar","kind"=>"call","operation"=>operation_record,
                    "args"=>[_encode_record(Num(a),w,path*"/arg") for a in Symbolics.arguments(raw)],"metadata"=>_scalar_metadata(raw,w,path))
    end
    if raw isa Number
        return Dict("tag"=>"scalar","kind"=>"constant","value"=>_encode_record(raw,w,path))
    end
    is_parameter(x) || _ioerror(path,"only atomic scalar variables and supported scalar ASTs are portable")
    record=Dict("tag"=>"scalar","kind"=>"variable","name"=>string(Symbolics.getname(raw)),"type"=>string(Symbolics.symtype(raw)),"metadata"=>_scalar_metadata(raw,w,path))
    key=(record["name"],record["type"]);canonical=_canonical(record)
    haskey(w.variables,key) && w.variables[key]!=canonical && _ioerror(path,"conflicting symbolic variable metadata for $(key[1])")
    w.variables[key]=canonical
    return record
end
function _decode_scalar(d,r,path)
    kind=d["kind"]
    if kind=="variable"
        key=(d["name"],d["type"]);canonical=_canonical(d)
        if haskey(r.variables,key)
            r.variables[key][1]==canonical || _ioerror(path,"conflicting symbolic variable declarations")
            return r.variables[key][2]
        end
        T=d["type"]=="Real" ? Real : d["type"]=="Complex" ? Complex : get(_IO_NUMBERS,d["type"],nothing)
        isnothing(T) && _ioerror(path,"unsupported variable scalar type")
        variable=Symbolics.variable(Symbol(d["name"]);T)
        raw=Symbolics.unwrap(variable)
        for (id,value) in get(d,"metadata",Dict())
            if startswith(id,"custom:")
                key_type=_decode_record(value["key"],r,path)
                key_type isa DataType || _ioerror(path,"metadata codec must restore a DataType key")
                raw=SymbolicUtils.setmetadata(raw,key_type,_decode_record(value["value"],r,path))
            else
                haskey(_IO_VARIABLE_METADATA,id) || _ioerror(path,"unknown variable metadata")
                raw=SymbolicUtils.setmetadata(raw,_IO_VARIABLE_METADATA[id],_decode_record(value,r,path))
            end
        end
        variable=Num(raw)
        r.variables[key]=(canonical,variable)
        return variable
    end
    kind=="constant" && return Num(_decode_record(d["value"],r,path))
    kind=="call" || _ioerror(path,"unknown symbolic node")
    operation=_decode_record(d["operation"],r,path*"/operation")
    args=[_decode_record(a,r,path) for a in d["args"]]
    raw=Symbolics.maketerm(typeof(Symbolics.unwrap(Symbolics.variable(:_io_scalar))),operation,Symbolics.unwrap.(args),nothing)
    for (id,value) in get(d,"metadata",Dict())
        if startswith(id,"custom:")
            key_type=_decode_record(value["key"],r,path)
            key_type isa DataType || _ioerror(path,"metadata codec must restore a DataType key")
            raw=SymbolicUtils.setmetadata(raw,key_type,_decode_record(value["value"],r,path))
        else
            haskey(_IO_VARIABLE_METADATA,id) || _ioerror(path,"unknown expression metadata")
            raw=SymbolicUtils.setmetadata(raw,_IO_VARIABLE_METADATA[id],_decode_record(value,r,path))
        end
    end
    return Num(raw)
end
