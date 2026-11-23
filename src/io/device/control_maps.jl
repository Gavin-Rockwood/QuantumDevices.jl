# Builtins save immutable definitions only, never compiled closures or preparation callables.
_control_map_record(map::AbstractControlMap,w,path) = _custom_record(map,w,path)
function _control_map_record(map::AffineControlMap,w,path)
    fields = (;name=map.name,inputs=map.inputs,outputs=map.outputs,
        response=reduce(vcat,(permutedims(collect(row)) for row in map.response)),offset=map.offset,
        input_bounds=map.input_bounds,output_bounds=map.output_bounds,metadata=map.metadata)
    return Dict("tag"=>"control_map","kind"=>"affine","version"=>1,"fields"=>_encode_record(fields,w,path))
end
function _control_map_record(map::TabulatedControlMap,w,path)
    fields = (;name=map.name,input=only(map.inputs),output=only(map.outputs),
        voltage_knots=map.voltage_knots,flux_knots=map.flux_knots,
        derivation=map.derivation,metadata=map.metadata,extrapolation=map.extrapolation)
    return Dict("tag"=>"control_map","kind"=>"tabulated","version"=>1,"fields"=>_encode_record(fields,w,path))
end
function _restore_control_map(d,r,path)
    get(d,"version",nothing)==1 || _ioerror(path,"unsupported control map version")
    f = _decode_record(d["fields"],r,path)
    f isa NamedTuple || _ioerror(path,"invalid control map fields")
    if d["kind"]=="affine"
        Set(keys(f))==Set((:name,:inputs,:outputs,:response,:offset,:input_bounds,:output_bounds,:metadata)) || _ioerror(path,"affine map field schema mismatch")
        return AffineControlMap(f.name;inputs=f.inputs,outputs=f.outputs,response=f.response,
            offset=f.offset,input_bounds=f.input_bounds,output_bounds=f.output_bounds,metadata=f.metadata)
    elseif d["kind"]=="tabulated"
        Set(keys(f))==Set((:name,:input,:output,:voltage_knots,:flux_knots,:derivation,:metadata,:extrapolation)) || _ioerror(path,"tabulated map field schema mismatch")
        return TabulatedControlMap(f.name;input=f.input,output=f.output,voltage_knots=f.voltage_knots,
            flux_knots=f.flux_knots,derivation=f.derivation,metadata=f.metadata,extrapolation=f.extrapolation)
    end
    _ioerror(path,"unknown control map kind")
end
# Schema 1 had direct controls only. Add empty fields without consulting registries.
function _migrate_control_map_records!(record)
    if record isa AbstractDict
        if get(record,"kind",nothing) in ("Leaf","Combine","CompositeRecipe") && haskey(record,"fields") && !haskey(record["fields"],"control_maps")
            record["fields"]["control_maps"] = Dict("tag"=>"dict","entries"=>Any[])
        end
        foreach(_migrate_control_map_records!,values(record))
    elseif record isa AbstractVector
        foreach(_migrate_control_map_records!,record)
    end
    return record
end
