# Enumerate only supported physical sources, never arbitrary metadata or solver state.
function parameter_inventory_roots(device)
    roots = NamedTuple[]
    function addroot(source,path,recipe;context=device,reference=nothing,parameters=nothing,terms=())
        components = isnothing(context) ? Dict() : context.components
        recipes = isnothing(context) ? Dict() : context.recipes
        push!(roots,(;source,path,recipe,components,recipes,reference,parameters,terms,context))
    end
    function calibration(result,path)
        prepared = result.provenance.prepared
        addroot(:calibration,path,result.recipe.model_recipe;context=prepared.context,
            reference=prepared.model.recipe,parameters=prepared.parameters,terms=result.recipe.terms)
    end
    function gates(repository,path,active)
        repository in active && throw(ArgumentError("gate repository cycle at $path"))
        push!(active,repository)
        try
            for (name,value) in sort!(collect(repository);by=p->string(first(p)))
                address = (path...,name)
                if value isa AbstractDict
                    gates(value,address,active)
                    continue
                end
                value isa QuantumGate || throw(ArgumentError("invalid gate at $address"))
                addroot(:gate,address,value.recipe.model_recipe;terms=value.recipe.terms)
                isnothing(value.calibration) || calibration(value.calibration,address)
                for (key,artifact) in sort!(collect(value.artifacts);by=p->string(first(p)))
                    location = (address...,:artifacts,key)
                    if artifact isa CalibrationResult
                        calibration(artifact,location)
                    elseif artifact isa NumericalModel
                        addroot(:model,location,artifact.recipe.definition;context=nothing,reference=artifact.recipe)
                    elseif artifact isa GateEvolution
                        reference = artifact.reference
                        addroot(:evolution,location,artifact.recipe.model_recipe;context=nothing,
                            reference=reference.recipe,parameters=get(reference,:parameters,Dict()),terms=artifact.recipe.terms)
                    end
                end
            end
        finally
            delete!(active,repository)
        end
    end
    for (name,recipe) in sort!(collect(device.recipes);by=p->string(first(p)))
        addroot(:recipe,(name,),recipe)
    end
    gates(device.gates,(),Base.IdSet())

    # Derivations own a map-free physical context and an explicit held input map.
    function derivation(map,address,active)
        map isa TabulatedControlMap || return
        data = map.derivation
        data isa NamedTuple && all(k->haskey(data,k),(:recipe,:context,:parameters,:parameter)) || return
        map in active && throw(ArgumentError("control-map derivation cycle at $address"))
        push!(active,map)
        try
            addroot(:control_map,address,data.recipe;context=data.context,parameters=data.parameters)
            maps_in_recipe(data.recipe,data.context.recipes,address,Any[],active)
        finally
            delete!(active,map)
        end
    end
    function maps_in_recipe(source,recipes,address,ancestry,map_active)
        source = try
            _parameter_definition(source,recipes)
        catch error
            error isa ArgumentError || rethrow()
            throw(ArgumentError("parameter root $address: $(error.msg)"))
        end
        source isa ResolvedRecipe && (source=source.definition)
        any(x->x===source,ancestry) && throw(ArgumentError("recipe cycle at $address"))
        ancestry = [ancestry;source]
        for (name,map) in sort!(collect(source.control_maps);by=p->string(first(p)))
            derivation(map,(address...,:map,name),map_active)
        end
        for input in recipe_inputs(source)
            maps_in_recipe(input.recipe,recipes,(address...,:input,input.name),ancestry,map_active)
        end
    end
    for (name,map) in sort!(collect(device.control_maps);by=p->string(first(p)))
        derivation(map,(:device,name),Base.IdSet())
    end
    for root in copy(roots)
        source = isnothing(root.reference) ? root.recipe : root.reference.definition
        maps_in_recipe(source,root.recipes,(root.source,root.path...),Any[],Base.IdSet())
    end
    return roots
end

# Discovery remains usable when effective values are missing or contradictory.
function collect_inventory_variables(device,roots)
    variables = Dict{String,Num}()
    origins = Dict{String,Any}()
    function include_inputs(inputs,path)
        for variable in inputs
            is_parameter(variable) || throw(ArgumentError("non-atomic input at $path"))
            name = string(variable)
            haskey(variables,name) && !isequal(variables[name],variable) &&
                throw(ArgumentError("distinct input symbols named $name at $(origins[name]) and $path"))
            variables[name] = variable
            get!(origins,name,path)
        end
    end
    function declarations(context,path)
        isnothing(context) && return
        for component in values(context.components); include_inputs(values(component.variables),path); end
        for term in values(context.interactions)
            validate_interaction_inputs(term)
            include_inputs(values(term.variables),path)
        end
    end
    declarations(device,(:device,))
    for root in roots
        declarations(root.context,(root.source,root.path...))
        inputs = try
            recipe_input_variables(root.recipe,root.components;extra_terms=root.terms,
                reference=root.reference,recipes=root.recipes)
        catch error
            error isa ArgumentError || rethrow()
            throw(ArgumentError("$(root.source) $(root.path): $(error.msg)"))
        end
        include_inputs(inputs,(root.source,root.path...))
    end
    return variables
end

"""Discover full symbolic input names across declarations and supported recipe contexts."""
device_variables(device) = collect_inventory_variables(device,parameter_inventory_roots(device))

# Read snapshots are independent; dictionary equality groups numeric keys with isequal.
function device_parameter_view(device)
    roots = parameter_inventory_roots(device)
    variables = collect_inventory_variables(device,roots)
    inventory = Dict(name=>Dict{Union{Missing,Number},Vector{NamedTuple}}() for name in keys(variables))
    validate_component_input_defaults(device.components)
    for root in roots
        overrides = root.source===:control_map && !isnothing(root.parameters) ? root.parameters : Dict()
        expressions = [physical_expression(term.expression) for term in root.terms]
        resolution = try
            resolve_recipe_inputs(root.recipe,root.components;overrides,extra_terms=root.terms,
                expressions,reference=root.reference,reference_parameters=root.parameters,recipes=root.recipes)
        catch error
            error isa ArgumentError || rethrow()
            throw(ArgumentError("$(root.source) $(root.path): $(error.msg)"))
        end
        gate_root = root.source in (:gate,:calibration,:evolution)
        for entry in resolution.nodes
            path = gate_root ? (:model,entry.path...) : entry.path
            append_input_locations!(inventory,entry.values,(source=root.source,path=root.path,node=path))
        end
        !gate_root || append_input_locations!(inventory,resolution.values,(source=root.source,path=root.path,node=()))
    end
    return inventory
end
function append_input_locations!(inventory,values,location)
    for (variable,value) in values
        locations = get!(inventory[string(variable)],value,NamedTuple[])
        location in locations || push!(locations,location)
    end
    return nothing
end
