# Old default repositories are interpreted only while upgrading verified records.
function _migrate_parameter_records!(definition,recipes,gates;schema,report)
    schema <= 3 || return nothing
    originals = copy(recipes)
    reader = RecordReader(recipes,Dict(),Dict(),Set(),Dict(),Dict())
    writer = RecordWriter()
    decode(record,path) = _decode_record(record,reader,path)
    encode(value,path) = _encode_record(value,writer,path)
    entries(record) = _legacy_entries(record,"parameter migration")
    items(record) = get(record,"values",Any[])
    pending = NamedTuple[]
    active = Set{String}()
    contexts = Dict{String,Any}()
    function term_defaults(term,path,used)
        fields = _legacy_fields(term)
        result = Dict{Num,Number}()
        if !isnothing(fields) && haskey(fields,"parameters")
            variables = decode(fields["variables"],path)
            for (name,value) in decode(fields["parameters"],path)
                variable = variables[name]
                result[variable] = value
                push!(pending,(;path,variable,value,used))
            end
        end
        return result
    end
    function context_values(record,path)
        fields = _legacy_fields(record)
        isnothing(fields) && return nothing
        all(key->haskey(fields,key),("name","components","interactions","recipes","metadata")) || return nothing
        components = _legacy_components(record,path)
        isnothing(components) && return nothing
        fingerprint = _canonical(record)
        haskey(contexts,fingerprint) && return contexts[fingerprint]
        used = Set{Num}()
        defaults = Dict{Num,Number}()
        groups = Any[]
        if haskey(fields,"explicit_parameters")
            values = Dict{Num,Number}(decode(fields["explicit_parameters"],path))
            push!(groups,("$path/device defaults",values));merge!(defaults,values)
            for (variable,value) in values; push!(pending,(;path,variable,value,used)); end
        end
        for component in Base.values(components)
            c = _legacy_fields(component)
            variables = decode(c["variables"],path); values = decode(c["parameters"],path)
            declared = Dict{Num,Number}(variables[name]=>value for (name,value) in values)
            push!(groups,("$path/component",declared));merge!(defaults,declared)
        end
        if haskey(fields,"interactions")
            for (_,term) in entries(fields["interactions"])
                values = term_defaults(term,path,used)
                push!(groups,("$path/interaction",values));merge!(defaults,values)
            end
        end
        validate_default_sources(groups)
        return get!(contexts,fingerprint,(;defaults,components,used,recipes=fields["recipes"]))
    end
    function named_source(source,context,path)
        if get(source,"tag",nothing)=="symbol"
            for (key,value) in entries(context.recipes)
                key["value"]==source["value"] && return value
            end
            _ioerror(path,"unknown named recipe $(source["value"])")
        end
        return source
    end
    function recipe_record(source)
        if get(source,"tag",nothing)=="recipe_ref"
            id=source["id"]
            haskey(originals,id) || _ioerror("parameter migration","unknown recipe node $id")
            return originals[id]
        end
        return source
    end
    function collect_defaults(source,context,path,seen=Set{String}())
        source=named_source(source,context,path)
        if get(source,"tag",nothing)=="recipe_ref"
            id=source["id"];id in seen && _ioerror(path,"recipe cycle")
            seen=union(seen,Set([id]))
        end
        record=recipe_record(source); fields=record["fields"]
        groups=Any[]
        if haskey(fields,"interactions")
            for term in items(fields["interactions"])
                push!(groups,(path,term_defaults(term,path,context.used)))
            end
        end
        for child in raw_inputs(record)
            append!(groups,collect_defaults(child,context,path,seen))
        end
        return groups
    end
    function raw_inputs(record)
        fields=record["fields"]
        record["kind"]=="Leaf" && return Any[]
        inputs=fields["inputs"]
        bindings=get(inputs,"tag",nothing)=="named_tuple" ? inputs["values"] : items(inputs)
        return [binding["fields"]["recipe"] for binding in bindings]
    end
    function root_symbols(source,context,path;reference=nothing,terms=())
        known=Set{Num}(); complete=true; ancestry=Set{String}()
        function saved_inputs(saved)
            fields=_legacy_fields(saved)
            union!(known,keys(decode(fields["parameters"],path)))
            for input in items(fields["inputs"]);saved_inputs(input["value"]);end
        end
        function visit(input)
            input=named_source(input,context,path)
            id=get(input,"id",nothing)
            id in ancestry && _ioerror(path,"recipe cycle")
            isnothing(id) || push!(ancestry,id)
            record=recipe_record(input);fields=record["fields"]
            if record["kind"]=="Leaf"
                owner=fields["component"]["value"]
                if haskey(context.components,owner)
                    union!(known,Base.values(decode(_legacy_fields(context.components[owner])["variables"],path)))
                else
                    complete=false
                end
            else
                for term in items(fields["interactions"])
                    union!(known,Base.values(decode(term["fields"]["variables"],path)))
                end
                foreach(visit,raw_inputs(record))
            end
            isnothing(id) || delete!(ancestry,id)
        end
        isnothing(reference) ? visit(source) : saved_inputs(reference)
        for term in terms;union!(known,Base.values(decode(term["fields"]["variables"],path)));end
        return known,complete
    end
    function materialize(source,context,defaults,inherited,path;reference=nothing,known=nothing,complete=true)
        source=named_source(source,context,path)
        saved=isnothing(reference) ? nothing : _legacy_fields(reference)
        captured=isnothing(saved) ? Dict{Num,Number}() : Dict{Num,Number}(decode(saved["parameters"],path))
        upstream=isnothing(saved) ? Any[] : items(saved["inputs"])
        ref=get(source,"tag",nothing)=="recipe_ref"
        id=ref ? source["id"] : nothing
        ref && id in active && _ioerror(path,"recipe cycle")
        ref && push!(active,id)
        try
            record=deepcopy(recipe_record(source));fields=record["fields"]
            own=Dict{Num,Number}(decode(fields["parameters"],path))
            if !isnothing(known)
                for (variable,value) in collect(own)
                    variable in known && continue
                    complete || _ioerror(path,"cannot prove old override $variable inactive without its historical context; assign inputs explicitly and resave")
                    delete!(own,variable)
                    push!(report,(kind=:dropped_default,path,variable=string(variable),value,reason=:inactive_override))
                end
            end
            assigned=merge(own,inherited)
            values=Dict{Num,Union{Missing,Number}}();origins=Dict{Num,Tuple}()
            local_inputs=Set{Num}()
            component_defaults=Dict{Num,Number}()
            if record["kind"]=="Leaf"
                owner=fields["component"]["value"]
                if haskey(context.components,owner)
                    c=_legacy_fields(context.components[owner])
                    variables=decode(c["variables"],path);idle=decode(c["parameters"],path)
                    union!(local_inputs,Base.values(variables))
                    component_defaults=Dict{Num,Number}(variables[name]=>value for (name,value) in idle)
                end
                isnothing(saved) || union!(local_inputs,keys(captured))
            else
                inputs=fields["inputs"]
                bindings=get(inputs,"tag",nothing)=="named_tuple" ? inputs["values"] : items(inputs)
                for binding in bindings
                    child= binding["fields"]["recipe"]
                    alias = binding["fields"]["name"]["value"]
                    reference_child = nothing
                    for entry in upstream
                        entry["key"]["value"]==alias && (reference_child=entry["value"])
                    end
                    updated,child_values=materialize(child,context,defaults,assigned,"$path -> $alias";reference=reference_child,known,complete)
                    binding["fields"]["recipe"]=updated
                    merge_input_values!(values,origins,child_values,(path,))
                end
                for term in items(fields["interactions"])
                    union!(local_inputs,Base.values(decode(term["fields"]["variables"],path)))
                end
            end
            for variable in local_inputs
                value=get(captured,variable,get(assigned,variable,get(values,variable,get(defaults,variable,Base.missing))))
                values[variable]=value
                if !haskey(assigned,variable) && !ismissing(value) && !haskey(origins,variable)
                    if !haskey(component_defaults,variable) || !isequal(component_defaults[variable],value)
                        own[variable]=value
                    end
                    push!(context.used,variable)
                end
            end
            isnothing(saved) || merge!(values,captured)
            fields["parameters"]=encode(own,path)
            fields["control_maps"]=upgrade(fields["control_maps"],context,"$path/control_maps")
            if ref
                new_id=recipe_fingerprint(record);recipes[new_id]=record
                updated=deepcopy(source);updated["id"]=new_id
                return updated,values
            end
            return record,values
        finally
            ref && delete!(active,id)
        end
    end
    empty_context=(;defaults=Dict{Num,Number}(),components=Dict(),used=Set{Num}(),recipes=Dict("tag"=>"dict","entries"=>Any[]))
    function upgrade(record,context,path;reference=nothing,parameters=nothing,captured=false)
        record isa AbstractVector && return Any[upgrade(v,context,path;reference,parameters,captured) for v in record]
        record isa AbstractDict || return record
        result=deepcopy(record)
        fields=_legacy_fields(record)
        if !isnothing(fields) && haskey(fields,"context")
            saved_context=context_values(fields["context"],path)
            if !isnothing(saved_context)
                context=saved_context
                captured=true
            end
        end
        local_context=context_values(record,path)
        if !isnothing(local_context)
            context=local_context
            reference=nothing
            parameters=nothing
        end
        tag=get(record,"tag",nothing)
        kind=get(record,"kind",nothing)
        if tag=="domain" && kind=="CalibrationResult"
            provenance=_legacy_fields(fields["provenance"])
            prepared_fields=_legacy_fields(provenance["prepared"])
            saved_context=context_values(prepared_fields["context"],path)
            if !isnothing(saved_context)
                context=saved_context;captured=true
            end
            reference=_legacy_fields(prepared_fields["model"])["recipe"]
            parameters=Dict{Num,Number}(decode(prepared_fields["parameters"],path))
        elseif tag=="domain" && kind=="PreparedGate"
            reference=_legacy_fields(fields["model"])["recipe"]
            parameters=Dict{Num,Number}(decode(fields["parameters"],path))
        elseif tag=="domain" && kind in ("NumericalModel","GateEvolution") && !captured
            context=empty_context
        end
        if tag=="domain" && kind=="GateEvolution"
            evolution_fields=_legacy_fields(fields["reference"])
            reference=evolution_fields["recipe"]
            stored=get(evolution_fields,"parameters",_legacy_fields(reference)["parameters"])
            parameters=Dict{Num,Number}(decode(stored,path))
        end
        if tag=="domain" && kind=="ResolvedRecipe"
            reference=record
        end
        if tag=="recipe_ref"
            haskey(originals,record["id"]) || return result
            groups=collect_defaults(record,context,path)
            validate_default_sources(vcat([(path,context.defaults)],groups))
            defaults=merge(context.defaults,Dict{Num,Number}(v=>x for (_,group) in groups for (v,x) in group))
            known,complete=root_symbols(record,context,path;reference)
            result,_=materialize(record,context,defaults,Dict(),path;reference,known,complete)
        elseif tag=="domain" && kind=="GateRecipe"
            f=result["fields"]
            groups=collect_defaults(f["model_recipe"],context,path)
            for term in items(f["terms"]);push!(groups,(path,term_defaults(term,path,context.used)));end
            validate_default_sources(vcat([(path,context.defaults)],groups))
            defaults=merge(context.defaults,Dict{Num,Number}(v=>x for (_,group) in groups for (v,x) in group))
            known,complete=root_symbols(f["model_recipe"],context,path;reference,terms=items(f["terms"]))
            endpoint,model_values=materialize(f["model_recipe"],context,defaults,Dict(),path;reference,known,complete)
            inserted=false
            # Gate-only removed defaults belong to the gate's own endpoint copy.
            current=recipes[endpoint["id"]]
            assignment=Dict{Num,Number}(decode(current["fields"]["parameters"],path))
            for term in items(f["terms"]), variable in Base.values(decode(term["fields"]["variables"],path))
                if !haskey(model_values,variable) && !haskey(assignment,variable)
                    value = isnothing(parameters) ? get(defaults,variable,missing) : get(parameters,variable,missing)
                    if !ismissing(value)
                        assignment[variable]=value;push!(context.used,variable);inserted=true
                    end
                end
            end
            current=deepcopy(current);current["fields"]["parameters"]=encode(assignment,path)
            new_id=recipe_fingerprint(current);recipes[new_id]=current
            endpoint["id"]=new_id
            if get(f["model_recipe"],"tag",nothing)!="symbol" || inserted
                f["model_recipe"]=endpoint
            end
        end
        if tag=="domain" && kind=="GateEvolution"
            evolution_reference=result["fields"]["reference"]
            saved=_legacy_fields(evolution_reference)
            if !isnothing(saved) && !haskey(saved,"parameters")
                model_reference=_legacy_fields(saved["recipe"])
                saved["parameters"]=deepcopy(model_reference["parameters"])
                _set_legacy_fields!(evolution_reference,saved)
                gate_fields=_legacy_fields(result["fields"]["recipe"])
                for term in items(gate_fields["terms"]), variable in Base.values(decode(term["fields"]["variables"],path))
                    haskey(parameters,variable) || push!(report,(kind=:unavailable_baseline,path,variable=string(variable),value=missing,reason=:legacy_evolution_gate_inputs))
                end
            end
        end
        fields=_legacy_fields(result)
        if !isnothing(local_context) && haskey(fields,"explicit_parameters")
            delete!(fields,"explicit_parameters");_set_legacy_fields!(result,fields)
        elseif tag=="domain" && kind=="Interaction"
            delete!(result["fields"],"parameters")
        end
        if tag=="domain" && kind=="GateRecipe"
            result["fields"]=Dict(key=>(key=="model_recipe" ? value : upgrade(value,context,"$path/$key";reference,parameters,captured)) for (key,value) in result["fields"])
            return result
        end
        return Dict(k=>upgrade(v,context,"$path/$k";reference,parameters,captured) for (k,v) in result)
    end
    migrated=upgrade(definition,empty_context,"device")
    base=context_values(definition,"device")
    new_gates=Dict(id=>upgrade(record,isnothing(base) ? empty_context : base,"gates/$id") for (id,record) in gates)
    # Reader IDs are temporary lookup identities after structural stripping, not
    # claims of canonical hashes. A schema-4 save re-interns every current record.
    # Newly interned recipe records also need their structural Interaction field upgrade.
    function strip!(record)
        if record isa AbstractDict
            get(record,"tag",nothing)=="domain" && get(record,"kind",nothing)=="Interaction" && delete!(record["fields"],"parameters")
            foreach(strip!,Base.values(record))
        elseif record isa AbstractVector
            foreach(strip!,record)
        end
        return record
    end
    foreach(strip!,Base.values(recipes))
    empty!(definition);merge!(definition,migrated)
    empty!(gates);merge!(gates,new_gates)
    for entry in pending
        entry.variable in entry.used && continue
        diagnostic=(kind=:dropped_default,path=entry.path,variable=string(entry.variable),value=entry.value,reason=:unused_or_overridden)
        diagnostic in report || push!(report,diagnostic)
    end
    return nothing
end
