# Custom codecs name implementations explicitly; portable records never carry a Julia closure.
const _IO_CODECS = Dict{String,NamedTuple}()
const _IO_SOURCES = Dict{String,NamedTuple}()
const _IO_MODULES = Dict{Symbol,String}()
"""
    register_io_codec!(id, T; version, encode, decode, source_id=nothing)

Register an explicit portable configuration codec. Callbacks are trusted local
code; bundled source is executed only by an opted-in folder load. Unknown values
never fall back to arbitrary Julia serialization.
"""
function register_io_codec!(id::AbstractString, ::Type{T}; version::Integer, encode, decode, source_id=nothing) where T
    source_id isa Union{Nothing,AbstractString} || throw(ArgumentError("source_id must be a string or nothing"))
    isempty(id) && throw(ArgumentError("codec ID must be nonempty"))
    version > 0 || throw(ArgumentError("codec version must be positive"))
    haskey(_IO_CODECS,id) && throw(ArgumentError("codec already registered: $id"))
    any(c -> c.type === T,values(_IO_CODECS)) && throw(ArgumentError("type already has a codec: $T"))
    _IO_CODECS[String(id)] = (;type=T,version=Int(version),encode,decode,source_id=isnothing(source_id) ? nothing : String(source_id))
    return nothing
end
_ioerror(path,message)=throw(ArgumentError("$path: $message"))
function _custom_record(x,w,path)
    for (id,codec) in _IO_CODECS
        x isa codec.type || continue
        !isnothing(codec.source_id) && push!(w.sources,String(codec.source_id))
        record=Dict("tag"=>"custom","id"=>id,"version"=>codec.version,
            "config"=>_encode_record(Base.invokelatest(codec.encode,x),w,path*"/config"))
        if !isnothing(codec.source_id)
            source_id=String(codec.source_id)
            source=_writer_source!(w,source_id)
            record["source"]=Dict("id"=>source_id,"fingerprint"=>_source_fingerprint(source))
        end
        return record
    end
    _ioerror(path,"unsupported $(typeof(x)); register an explicit IO codec")
end
