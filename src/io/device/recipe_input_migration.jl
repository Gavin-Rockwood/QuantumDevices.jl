# Bind legacy leaf declarations per owning physical context after old IDs are checked.
function _migrate_recipe_input_records!(definition,recipes,gates;schema)
    schema < 5 || return nothing
    originals = copy(recipes)
    reader = RecordReader(recipes,Dict(),Dict(),Set(),Dict(),Dict())
    writer = RecordWriter()
    active = Set{String}()
    function upgrade(record,components,path;historical=false,held=nothing)
        record isa AbstractVector && return Any[upgrade(v,components,path;historical,held) for v in record]
        record isa AbstractDict || return record
        fields = _legacy_fields(record)
        if !isnothing(fields) && haskey(fields,"context")
            own = _legacy_components(fields["context"],path)
            isnothing(own) || (components=own; historical=true)
        end
        own = _legacy_components(record,path)
        isnothing(own) || (components=own)
        kind = get(record,"kind",nothing)
        if get(record,"tag",nothing)=="domain" && kind in ("NumericalModel","GateEvolution") && !historical
            components=nothing
        end
        if get(record,"tag",nothing)=="domain" && kind=="ResolvedRecipe"
            held = _decode_record(fields["parameters"],reader,path)
        end
        result = deepcopy(record)
        if get(record,"tag",nothing)=="recipe_ref"
            id=record["id"]
            id in active && _ioerror(path,"cyclic recipe inputs")
            haskey(originals,id) || _ioerror(path,"unknown recipe node")
            push!(active,id)
            migrated = try
                upgrade(originals[id],components,"$path/recipe";historical,held)
            finally
                delete!(active,id)
            end
            new_id=recipe_fingerprint(migrated)
            recipes[new_id]=migrated
            result["id"]=new_id
            return result
        end
        if kind=="Leaf" && haskey(record,"fields") && !haskey(record["fields"],"input_defaults")
            leaf=record["fields"]
            name=leaf["component"]["value"]
            defaults=nothing
            if !isnothing(components) && haskey(components,name)
                component=_legacy_fields(components[name])
                variables=_decode_record(component["variables"],reader,path)
                idle=_decode_record(component["parameters"],reader,path)
                defaults=RecipeMap(Dict{Num,Union{Missing,Number}}(v=>get(idle,k,missing) for (k,v) in variables))
            elseif !isnothing(held)
                defaults=RecipeMap(Dict{Num,Union{Missing,Number}}(v=>missing for v in keys(held)))
            end
            result["fields"]["input_defaults"]=_encode_record(defaults,writer,path)
        end
        return Dict(k=>upgrade(v,components,"$path/$k";historical,held) for (k,v) in result)
    end
    updated=upgrade(definition,nothing,"device")
    components=_legacy_components(definition,"device")
    updated_gates=Dict(id=>upgrade(record,components,"gates/$id") for (id,record) in gates)
    empty!(definition);merge!(definition,updated)
    empty!(gates);merge!(gates,updated_gates)
    return nothing
end
