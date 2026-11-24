# Upgrade legacy declaration records before constructing current runtime objects.
# All original recipe content IDs have already been verified by the folder loader.
function _legacy_fields(record)
    record isa AbstractDict || return nothing
    tag = get(record,"tag",nothing)
    tag == "domain" && return record["fields"]
    tag == "named_tuple" && return Dict(zip(record["keys"],record["values"]))
    return nothing
end
function _set_legacy_fields!(record, fields)
    if get(record,"tag",nothing) == "named_tuple"
        names = [name for name in record["keys"] if haskey(fields,name)]
        append!(names,sort!([name for name in keys(fields) if !(name in names)]))
        record["keys"] = names
        record["values"] = Any[fields[name] for name in names]
    else
        record["fields"] = fields
    end
    return record
end
function _legacy_entries(record, path)
    record isa AbstractDict && get(record,"tag",nothing) in ("dict","recipe_map") ||
        _ioerror(path,"legacy declaration map is not a dictionary")
    return record["entries"]
end
function _legacy_generator_fields(record, path)
    get(record,"tag",nothing) == "domain" && get(record,"kind",nothing) == "OperatorGenerator" ||
        _ioerror(path,"legacy generator requires explicit function/default records; resave with a one-argument function")
    fields = record["fields"]
    Set(keys(fields)) == Set(("f","default_dimension")) || _ioerror(path,"invalid legacy generator fields")
    return fields
end
function _legacy_generator_dimension(record, path)
    value = _legacy_generator_fields(record,path)["default_dimension"]
    get(value,"tag",nothing) == "number" || _ioerror(path,"invalid legacy generator dimension")
    dimension = tryparse(Int,value["value"])
    isnothing(dimension) && _ioerror(path,"invalid legacy generator dimension")
    return validate_operator_dimension(dimension)
end
function _legacy_component_dimension(generators, path)
    dimensions = unique(_legacy_generator_dimension(entry[2],path) for entry in _legacy_entries(generators,path))
    return length(dimensions) == 1 ? only(dimensions) : nothing
end
_legacy_dimension_record(dimension) = isnothing(dimension) ? Dict("tag"=>"nothing") :
    Dict("tag"=>"number","type"=>"Int64","value"=>string(dimension))

# Context records can be live Component records or captured named tuples.
function _legacy_components(context, path)
    fields = _legacy_fields(context)
    isnothing(fields) && return nothing
    haskey(fields,"components") || return nothing
    map = fields["components"]
    get(map,"tag",nothing) in ("dict","recipe_map") || return nothing
    return Dict(entry[1]["value"] => entry[2] for entry in _legacy_entries(map,path))
end

function _migrate_leaf_generators!(record, components, path)
    fields = record["fields"]
    owner = fields["component"]["value"]
    component = isnothing(components) ? nothing : get(components,owner,nothing)
    # Keep legacy per-operator defaults until the effective old dimension is known.
    merged = Dict{String,Any}()
    if !isnothing(component)
        cfields = _legacy_fields(component)
        for (key,generator) in _legacy_entries(cfields["operator_generators"],path)
            merged[key["value"]] = generator
        end
    end
    for (key,generator) in _legacy_entries(fields["generators"],path)
        get(key,"tag",nothing) == "operator" && key["owner"] == owner || _ioerror(path,"invalid legacy generator override owner")
        merged[key["name"]] = generator
    end
    explicit = get(fields["dimension"],"tag",nothing) != "nothing"
    if !explicit
        isnothing(component) && _ioerror(path,"legacy implicit dimension needs component context; set Leaf.dimension and resave")
        dimensions = unique(_legacy_generator_dimension(generator,path) for generator in values(merged))
        length(dimensions) == 1 || _ioerror(path,"ambiguous legacy generator defaults; set Leaf.dimension and resave")
        fields["dimension"] = _legacy_dimension_record(only(dimensions))
    end
    if haskey(fields,"kwargs")
        kwargs = fields["kwargs"]
        get(kwargs,"tag",nothing) == "named_tuple" || _ioerror(path,"invalid legacy leaf kwargs")
        if !isempty(kwargs["keys"])
            isnothing(component) && _ioerror(path,"legacy keyword migration needs component context")
            all(generator -> get(_legacy_generator_fields(generator,path)["f"],"tag",nothing) == "builtin_generator", values(merged)) ||
                _ioerror(path,"custom legacy generator kwargs require an explicit one-argument closure; adapt and resave")
        end
        delete!(fields,"kwargs")
    end
    return record
end

function _migrate_operator_records!(definition, recipes, gate_records; schema)
    schema in (1,2) || return nothing
    original_recipes = copy(recipes)
    active = Set{String}()
    function upgrade(record, components, path; captured=false)
        record isa AbstractVector && return Any[upgrade(value,components,"$path/$i";captured) for (i,value) in enumerate(record)]
        record isa AbstractDict || return record
        result = deepcopy(record)
        # A nearby saved physical context determines recipe interpretation, regardless
        # of serialized field order. Never reuse the first device's defaults globally.
        fields = _legacy_fields(record)
        if !isnothing(fields) && haskey(fields,"context")
            local_components = _legacy_components(fields["context"],path)
            if !isnothing(local_components)
                components = local_components
                captured = true
            end
        end
        local_components = _legacy_components(record,path)
        isnothing(local_components) || (components = local_components)
        tag = get(record,"tag",nothing)
        # Historical numerical artifacts need their own captured physical context;
        # live device defaults may have changed since those arrays were built.
        if tag == "domain" && get(record,"kind",nothing) in ("NumericalModel","GateEvolution") && !captured
            components = nothing
        end
        if tag == "recipe_ref"
            id = record["id"]
            haskey(original_recipes,id) || _ioerror(path,"unknown legacy recipe node")
            id in active && _ioerror(path,"cyclic legacy recipe nodes")
            push!(active,id)
            upgraded = try
                upgrade(original_recipes[id],components,"$path/recipe";captured)
            finally
                delete!(active,id)
            end
            # Per-occurrence upgrades can differ under two physical contexts.
            new_id = recipe_fingerprint(upgraded)
            recipes[new_id] = upgraded
            result["id"] = new_id
            return result
        elseif tag == "builtin_generator"
            get(record,"version",1) == 1 || _ioerror(path,"unsupported legacy builtin generator version")
            id = get(record,"kind",nothing)
            haskey(BUILTIN_OPERATOR_FUNCTIONS,id) || _ioerror(path,"unknown legacy builtin generator")
            return Dict("tag"=>"operator_function","version"=>1,"id"=>id)
        elseif tag == "domain" && get(record,"kind",nothing) == "OperatorGenerator"
            return upgrade(_legacy_generator_fields(record,path)["f"],components,"$path/f";captured)
        end
        if get(record,"kind",nothing) == "Leaf" && haskey(record,"fields")
            _migrate_leaf_generators!(result,components,path)
        elseif !isnothing(fields) && haskey(fields,"operator_generators") && haskey(fields,"hamiltonian")
            dimension = _legacy_component_dimension(fields["operator_generators"],path)
            _set_legacy_fields!(result,merge(fields,Dict("dimension"=>_legacy_dimension_record(dimension))))
        end
        return Dict(key=>upgrade(value,components,"$path/$key";captured) for (key,value) in result)
    end
    upgraded_definition = upgrade(definition,nothing,"device")
    # Gate records own independent captured contexts; the device is only fallback.
    components = _legacy_components(definition,"device")
    upgraded_gates = Dict(id=>upgrade(record,components,"gates/$id") for (id,record) in gate_records)
    empty!(definition); merge!(definition,upgraded_definition)
    empty!(gate_records); merge!(gate_records,upgraded_gates)
    return nothing
end
