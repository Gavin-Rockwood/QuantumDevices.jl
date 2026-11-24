# Recreate legacy records from current fixtures, retaining valid content hashes.
function legacy_generator_records!(record; override_dimension=2)
    QD = QuantumDevices
    if record isa AbstractDict
        fields = QD._legacy_fields(record)
        if !isnothing(fields) && haskey(fields,"operator_generators") && haskey(fields,"hamiltonian")
            dimension = pop!(fields,"dimension")
            for entry in fields["operator_generators"]["entries"]
                entry[2] = Dict("tag"=>"domain","kind"=>"OperatorGenerator",
                    "fields"=>Dict("f"=>entry[2],"default_dimension"=>deepcopy(dimension)))
            end
            QD._set_legacy_fields!(record,fields)
        elseif get(record,"kind",nothing) == "Leaf" && haskey(record,"fields")
            fields = record["fields"]
            fields["kwargs"] = Dict("tag"=>"named_tuple","keys"=>String[],"values"=>Any[])
            for entry in fields["generators"]["entries"]
                entry[2] = Dict("tag"=>"domain","kind"=>"OperatorGenerator",
                    "fields"=>Dict("f"=>entry[2],"default_dimension"=>QD._legacy_dimension_record(override_dimension)))
            end
        elseif get(record,"tag",nothing) == "operator_function"
            record["tag"] = "builtin_generator"
            record["kind"] = pop!(record,"id")
        end
        foreach(value->legacy_generator_records!(value;override_dimension),values(record))
    elseif record isa AbstractVector
        foreach(value->legacy_generator_records!(value;override_dimension),record)
    end
    return record
end
