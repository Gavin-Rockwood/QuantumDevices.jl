# Analysis results are numerical data; decoding never propagates, tracks or refits.
const _IO_ANALYSIS_TYPES=Dict("FloquetBasis"=>FloquetBasis,"FloquetSweep"=>FloquetSweep,
    "AvoidedCrossingFit"=>AvoidedCrossingFit,"ResonanceResult"=>ResonanceResult)
function _analysis_record(value,w,path)
    T=_IO_ANALYSIS_TYPES[string(nameof(typeof(value)))]
    T((getfield(value,k) for k in fieldnames(T))...) # Stored-data validation only.
    return Dict("tag"=>"analysis","version"=>1,"kind"=>string(nameof(typeof(value))),"fields"=>_fields_record(value,w,path))
end
function _restore_analysis(record,r,path)
    record["version"]==1 || _ioerror(path,"unsupported analysis version")
    T=get(_IO_ANALYSIS_TYPES,record["kind"],nothing)
    isnothing(T) && _ioerror(path,"unknown analysis type")
    Set(keys(record["fields"]))==Set(String.(fieldnames(T))) || _ioerror(path,"analysis field schema mismatch")
    values=(_decode_record(record["fields"][String(k)],r,"$path/$k") for k in fieldnames(T))
    return T(values...)
end
