using SparseArrays
const _IO_NUMBERS = Dict(string(T)=>T for T in (Int8,Int16,Int32,Int64,Int128,UInt8,UInt16,UInt32,UInt64,UInt128,Float16,Float32,Float64,BigInt,BigFloat))
const _IO_OPERATIONS = Dict(string(f)=>f for f in (+,-,*,/,^,conj,real,imag,abs,abs2,sqrt,exp,log,sin,cos,tan,sinpi,cospi,sinh,cosh,asin,acos,atan))
mutable struct RecordWriter
    recipes::Dict{String,Any}
    arrays::Dict{String,Any}
    scopes::Dict{Symbol,String}
    active::Base.IdSet{Any}
    sources::Set{String}
    variables::Dict{Tuple{String,String},String}
    native::Dict{String,Any}
    source_inputs::Dict{String,Any}
    source_data::Dict{String,Any}
end
RecordWriter(;code_sources=Dict())=RecordWriter(Dict(),Dict(),Dict(),Base.IdSet(),Set(),Dict(),Dict(),Dict{String,Any}(code_sources),Dict())
mutable struct RecordReader
    recipes::Dict{String,Any}
    arrays::Dict{String,Any}
    scopes::Dict{String,Symbol}
    active::Set{String}
    variables::Dict{Tuple{String,String},Tuple{String,Num}}
    native::Dict{String,Any}
    containers::Base.IdSet{Any}
end
RecordReader(recipes,arrays,scopes,active,variables,native)=RecordReader(recipes,arrays,scopes,active,variables,native,Base.IdSet())
_scope(w,scope)=get!(w.scopes,scope) do; "s$(length(w.scopes)+1)";end
_unScope(r,id)=get!(r.scopes,id) do; gensym(:loaded_scope);end
# Canonical bytes ignore dictionary insertion order, while lists remain ordered.
function _canonical(x)
    x isa AbstractDict && return "{"*join((repr(k)*":"*_canonical(x[k]) for k in sort!(collect(keys(x)))),",")*"}"
    x isa AbstractVector && return "["*join(_canonical.(x),",")*"]"
    return repr(x)
end
recipe_fingerprint(record)=bytes2hex(SHA.sha256(_canonical(record)))
function _little_endian_bytes(array)
    bytes=collect(reinterpret(UInt8,vec(array)))
    if Base.ENDIAN_BOM != 0x04030201
        width=eltype(array)<:Complex ? sizeof(eltype(array))÷2 : sizeof(eltype(array))
        for start in 1:width:length(bytes);reverse!(@view bytes[start:start+width-1]);end
    end
    return bytes
end
function _array_record(x,w,path)
    a=Array(x)
    if isconcretetype(eltype(a)) && isbitstype(eltype(a)) && eltype(a) <: Union{Bool,Int8,Int16,Int32,Int64,Int128,UInt8,UInt16,UInt32,UInt64,UInt128,Float16,Float32,Float64,ComplexF32,ComplexF64}
        bytes=_little_endian_bytes(a)
        descriptor=Dict("type"=>string(eltype(a)),"size"=>collect(size(a)),"bytes"=>bytes2hex(bytes))
        id=recipe_fingerprint(descriptor)
        w.arrays[id]=bytes
        return Dict("tag"=>"array_payload","id"=>id,"type"=>string(eltype(a)),"size"=>collect(size(a)))
    end
    x in w.active && _ioerror(path,"cyclic array")
    push!(w.active,x)
    try
        return Dict("tag"=>"array","size"=>collect(size(a)),"values"=>[_encode_record(v,w,"$path/$i") for (i,v) in enumerate(vec(a))])
    finally
        delete!(w.active,x)
    end
end
function _encode_record(x,w,path)
    x===nothing && return Dict("tag"=>"nothing")
    x===missing && return Dict("tag"=>"missing")
    x isa Num && return _encode_scalar(x,w,path)
    x isa Bool && return Dict("tag"=>"bool","value"=>x)
    x isa BigFloat && return Dict("tag"=>"bigfloat","precision"=>precision(x),"value"=>string(x))
    x isa Real && haskey(_IO_NUMBERS,string(typeof(x))) && return Dict("tag"=>"number","type"=>string(typeof(x)),"value"=>string(x))
    x isa Rational && return Dict("tag"=>"rational","n"=>_encode_record(numerator(x),w,path),"d"=>_encode_record(denominator(x),w,path))
    x isa Complex && return Dict("tag"=>"complex","real"=>_encode_record(real(x),w,path),"imag"=>_encode_record(imag(x),w,path))
    x isa Symbol && return Dict("tag"=>"symbol","value"=>String(x))
    x isa AbstractString && return Dict("tag"=>"string","value"=>String(x))
    x isa Recipe && return intern_recipe!(w,x,path)
    x isa QuantumObject && return _quantum_record(x,w,path)
    x isa Sym && return Dict("tag"=>"sym","value"=>_encode_record(x.expr,w,path))
    x isa OperatorRef && return Dict("tag"=>"operator","owner"=>String(x.owner),"name"=>String(x.name))
    x isa BoundOperatorRef && return Dict("tag"=>"bound_operator","scope"=>_scope(w,x.scope),"input"=>String(x.input),"physical"=>_encode_record(x.physical,w,path))
    x isa SparseVector && return Dict("tag"=>"sparse_vector","length"=>length(x),"indices"=>_encode_record(x.nzind,w,path),"values"=>_encode_record(x.nzval,w,path))
    x isa SparseMatrixCSC && return Dict("tag"=>"sparse","size"=>collect(size(x)),"colptr"=>_encode_record(x.colptr,w,path),"rowval"=>_encode_record(x.rowval,w,path),"nzval"=>_encode_record(x.nzval,w,path))
    x isa NamedTuple && return Dict("tag"=>"named_tuple","keys"=>String.(collect(keys(x))),"values"=>[_encode_record(v,w,"$path/$k") for (k,v) in pairs(x)])
    x isa Pair && return Dict("tag"=>"pair","key"=>_encode_record(first(x),w,path),"value"=>_encode_record(last(x),w,path))
    x isa Tuple && return Dict("tag"=>"tuple","values"=>[_encode_record(v,w,"$path/$i") for (i,v) in enumerate(x)])
    x isa AbstractDict && return _dict_record(x,w,path)
    x isa AbstractSet && return Dict("tag"=>"set","values"=>sort!([_encode_record(v,w,path) for v in x];by=_canonical))
    x isa TrackedData && return _encode_domain_record(x,w,path)
    x isa AbstractArray && return _array_record(x,w,path)
    haskey(_IO_OPERATIONS,string(x)) && any(f->x===f,values(_IO_OPERATIONS)) && return Dict("tag"=>"operation","id"=>string(x))
    x===cubic_dressing_trajectory && return Dict("tag"=>"dressing_trajectory")
    x===gate_infidelity && return Dict("tag"=>"objective","id"=>"gate_infidelity")
    x isa SciMLBase.ObjSense && return Dict("tag"=>"objective_sense","value"=>string(x))
    x isa SciMLBase.NoAD && return Dict("tag"=>"no_ad")
    x isa SciMLBase.ReturnCode.T && return Dict("tag"=>"retcode","value"=>string(x))
    x isa BuiltinGenerator && return Dict("tag"=>"builtin_generator","version"=>1,"kind"=>String(x.kind))
    return _encode_domain_record(x,w,path)
end
function _dict_record(x,w,path)
    x in w.active && _ioerror(path,"cyclic container")
    push!(w.active,x)
    try
        entries=[[_encode_record(k,w,"$path/key"),_encode_record(v,w,"$path/$k")] for (k,v) in x]
        sort!(entries;by=e->_canonical(e[1]))
        return Dict("tag"=>x isa RecipeMap ? "recipe_map" : "dict","entries"=>entries)
    finally
        delete!(w.active,x)
    end
end
function _typed_collect(values)
    isempty(values) && return Any[]
    T=reduce(typejoin,typeof.(values))
    return T[values...]
end
function _decode_record_inner(d,r,path)
    d isa AbstractDict && haskey(d,"tag") || _ioerror(path,"invalid record")
    tag=d["tag"]
    tag=="nothing" && return nothing
    tag=="missing" && return missing
    tag=="bool" && return Bool(d["value"])
    tag=="bigfloat" && return setprecision(BigFloat,Int(d["precision"])) do; parse(BigFloat,d["value"]);end
    tag=="number" && return parse(get(_IO_NUMBERS,d["type"]) do; _ioerror(path,"unknown numeric type");end,d["value"])
    tag=="complex" && return complex(_decode_record(d["real"],r,path),_decode_record(d["imag"],r,path))
    tag=="rational" && return _decode_record(d["n"],r,path)//_decode_record(d["d"],r,path)
    tag=="symbol" && return Symbol(d["value"])
    tag=="string" && return String(d["value"])
    tag=="scalar" && return _decode_scalar(d,r,path)
    tag=="recipe_ref" && return restore_recipe_reference(r,d,path)
    tag=="operator" && return OperatorRef(Symbol(d["owner"]),Symbol(d["name"]))
    tag=="bound_operator" && return BoundOperatorRef(_unScope(r,d["scope"]),Symbol(d["input"]),_decode_record(d["physical"],r,path))
    tag=="sym" && return Sym(_decode_record(d["value"],r,path))
    tag=="pair" && return _decode_record(d["key"],r,path)=>_decode_record(d["value"],r,path)
    if tag in ("tuple","named_tuple","set","array")
        values=[_decode_record(v,r,"$path/$i") for (i,v) in enumerate(d["values"])]
        tag=="tuple" && return Tuple(values)
        tag=="named_tuple" && return NamedTuple{Tuple(Symbol.(d["keys"]))}(Tuple(values))
        tag=="set" && return Set(values)
        return reshape(_typed_collect(values),Tuple(Int.(d["size"])))
    end
    if tag in ("dict","recipe_map")
        entries=[_decode_record(e[1],r,path)=>_decode_record(e[2],r,path) for e in d["entries"]]
        value=isempty(entries) ? Dict() : Dict(entries)
        return tag=="recipe_map" ? _io_owned_map(value) : value
    end
    if tag=="array_payload"
        haskey(r.arrays,d["id"]) || _ioerror(path,"missing numerical payload")
        bytes=r.arrays[d["id"]]
        types=merge(_IO_NUMBERS,Dict("Bool"=>Bool,"ComplexF32"=>ComplexF32,"ComplexF64"=>ComplexF64))
        T=get(types,d["type"],nothing)
        !isnothing(T) && isbitstype(T) || _ioerror(path,"unsupported payload type")
        dims=Int.(d["size"])
        all(>=(0),dims) && big(prod(big.(dims)))*sizeof(T)==length(bytes) || _ioerror(path,"payload shape mismatch")
        T===Bool && !all(b->b in (0x00,0x01),bytes) && _ioerror(path,"invalid Boolean payload")
        if Base.ENDIAN_BOM != 0x04030201
            bytes=copy(bytes);width=T<:Complex ? sizeof(T)÷2 : sizeof(T)
            for start in 1:width:length(bytes);reverse!(@view bytes[start:start+width-1]);end
        end
        return copy(reshape(reinterpret(T,bytes),Tuple(dims)))
    end
    if tag=="sparse_vector"
        n=Int(d["length"]);indices=_decode_record(d["indices"],r,path);values=_decode_record(d["values"],r,path)
        n>=0 && length(indices)==length(values) && issorted(indices) && allunique(indices) && all(i->1<=i<=n,indices) || _ioerror(path,"invalid sparse vector")
        return SparseVector(n,indices,values)
    end
    if tag=="sparse"
        dims=Int.(d["size"]);cp=_decode_record(d["colptr"],r,path);rv=_decode_record(d["rowval"],r,path);nz=_decode_record(d["nzval"],r,path)
        length(cp)==dims[2]+1 && first(cp)==1 && last(cp)==length(nz)+1 && issorted(cp) && length(rv)==length(nz) && all(i->1<=i<=dims[1],rv) || _ioerror(path,"invalid sparse storage")
        return SparseMatrixCSC(dims...,cp,rv,nz)
    end
    tag=="quantum" && return _restore_quantum(d,r,path)
    tag=="operation" && return get(_IO_OPERATIONS,d["id"]) do; _ioerror(path,"unknown operation");end
    tag=="objective" && d["id"]=="gate_infidelity" && return gate_infidelity
    tag=="dressing_trajectory" && return cubic_dressing_trajectory
    if tag=="objective_sense"
        d["value"] in ("MinSense","MaxSense") || _ioerror(path,"unknown objective sense")
        return d["value"]=="MinSense" ? SciMLBase.MinSense : SciMLBase.MaxSense
    end
    tag=="no_ad" && return SciMLBase.NoAD()
    if tag=="retcode"
        name=Symbol(last(split(d["value"],'.')))
        name in propertynames(SciMLBase.ReturnCode) || _ioerror(path,"unknown return code")
        value=getproperty(SciMLBase.ReturnCode,name)
        value isa SciMLBase.ReturnCode.T || _ioerror(path,"not a return-code value")
        return value
    end
    if tag=="builtin_generator"
        get(d,"version",1)==1 || _ioerror(path,"unsupported builtin generator version")
        kind=Symbol(d["kind"])
        kind in (:pauli_x,:pauli_y,:pauli_z,:destroy,:number,:charge,:jump) || _ioerror(path,"unknown builtin generator")
        return BuiltinGenerator(kind)
    end
    if tag=="custom"
        codec=get(_IO_CODECS,d["id"],nothing)
        !isnothing(codec) && codec.version==d["version"] || _ioerror(path,"missing codec/version $(d["id"])")
        if haskey(d,"source")
            source=d["source"];id=source["id"]
            codec.source_id==id && haskey(_IO_SOURCES,id) && _source_fingerprint(_IO_SOURCES[id])==source["fingerprint"] || _ioerror(path,"custom source identity mismatch")
        end
        value=Base.invokelatest(codec.decode,_decode_record(d["config"],r,path*"/config"))
        value isa codec.type || _ioerror(path,"codec returned wrong type")
        return value
    end
    return _restore_domain(d,r,path)
end

# These values were decoded into an independent owner already. Freeze container
# structure without mistaking a registered mutable config for shared source data.
_io_owned(value::AbstractDict)=_io_owned_map(value)
_io_owned(value::AbstractVector)=_io_owned_vector(value)
_io_owned(value::Tuple)=map(_io_owned,value)
_io_owned(value::NamedTuple)=map(_io_owned,value)
_io_owned(value::Pair)=_io_owned(first(value))=>_io_owned(last(value))
_io_owned(value)=value
function _io_owned_map(source)
    items=Tuple(_io_owned(k)=>_io_owned(v) for (k,v) in source)
    K=isempty(items) ? keytype(source) : reduce(typejoin,typeof(first(p)) for p in items)
    V=isempty(items) ? valtype(source) : reduce(typejoin,typeof(last(p)) for p in items)
    return RecipeMap{K,V,typeof(items)}(items)
end
function _io_owned_vector(source)
    items=Tuple(_io_owned(v) for v in source)
    T=isempty(items) ? eltype(source) : reduce(typejoin,typeof.(items))
    return RecipeVector{T,typeof(items)}(items)
end

function _decode_record(d,r,path)
    d isa AbstractDict || _ioerror(path,"invalid record")
    d in r.containers && _ioerror(path,"cyclic logical record")
    push!(r.containers,d)
    try
        return _decode_record_inner(d,r,path)
    finally
        delete!(r.containers,d)
    end
end

function _encode_domain_record(x,w,path)
    x in w.active && _ioerror(path,"cyclic domain object")
    push!(w.active,x)
    try
        return _domain_record(x,w,path)
    finally
        delete!(w.active,x)
    end
end
