# Resolve physical membership without generating any component matrices.
function _recipe_map_variables(device, owners)
    return Num[v for owner in owners for v in values(lookup(device.components,owner,"component").variables)]
end
function _copy_recipe_maps(device, maps, owners, frozen)
    result = Dict{Symbol,AbstractControlMap}()
    for (name,map) in maps
        map isa AbstractControlMap && name === map.name || throw(ArgumentError("invalid recipe control map entry $name"))
        isnothing(device) || _validate_control_map_outputs(map,_recipe_map_variables(device,owners))
        result[name] = _copy_control_map(map)
    end
    return frozen ? RecipeMap(result) : result
end
# Local overrides settle inherited conflicts before automatic registry additions.
function _capture_recipe_maps(device, bindings, selection, owners)
    selection === :auto || !(selection isa Symbol) || throw(ArgumentError("control_maps must be :auto or a collection"))
    localmaps = Dict{Symbol,AbstractControlMap}()
    if selection !== :auto
        for item in selection
            map = if item isa Symbol
                isnothing(device) && throw(ArgumentError("named control maps require a device"))
                lookup(device.control_maps,item,"control map")
            else
                item
            end
            map isa AbstractControlMap || throw(ArgumentError("control_maps must contain maps or registered names"))
            haskey(localmaps,map.name) && throw(ArgumentError("duplicate local control map $(map.name)"))
            localmaps[map.name] = map
        end
    end
    result = Dict{Symbol,AbstractControlMap}()
    for binding in bindings
        source = binding.recipe isa Symbol ? lookup(device.recipes,binding.recipe,"recipe") : binding.recipe
        for (name,map) in source.control_maps
            haskey(localmaps,name) && continue
            haskey(result,name) && !isequal(result[name],map) && throw(ArgumentError("conflicting inherited control map $name; supply a local override"))
            result[name] = map
        end
    end
    merge!(result,localmaps)
    if selection === :auto && !isnothing(device)
        variables = _recipe_map_variables(device,owners)
        for name in sort!(collect(keys(device.control_maps)))
            map = device.control_maps[name]
            haskey(result,name) || all(out -> any(isequal(out),variables),map.outputs) && (result[name]=map)
        end
    end
    return _copy_recipe_maps(device,result,owners,true)
end

# Derivation recipe scopes are local binders, not calibration identity. Preserve
# physical fields, input aliases, maps and every derivation/context datum.
_control_map_recipe_value(value) = value
_control_map_recipe_value(value::Recipe) = _control_map_recipe_content(value)
_control_map_recipe_value(value::NamedTuple) = map(_control_map_recipe_value,value)
_control_map_recipe_value(value::Tuple) = map(_control_map_recipe_value,value)
_control_map_recipe_value(value::Pair) = _control_map_recipe_value(first(value)) => _control_map_recipe_value(last(value))
_control_map_recipe_value(value::AbstractDict) = Dict(_control_map_recipe_value(k)=>_control_map_recipe_value(v) for (k,v) in value)
_control_map_recipe_value(value::AbstractArray) = (;dimensions=size(value),items=Tuple(_control_map_recipe_value(v) for v in value))
_control_map_recipe_value(value::BoundOperatorRef) = (;kind=:bound_operator,input=value.input,physical=value.physical)
_control_map_recipe_value(value::Sym) = (;kind=:expression,content=_control_map_recipe_value(value.expr))
_control_map_recipe_value(value::InputBinding) = (;kind=:binding,name=value.name,
    recipe=_control_map_recipe_value(value.recipe),owners=_control_map_recipe_value(value.owners),
    operator_names=_control_map_recipe_value(value.operator_names))
_control_map_recipe_value(value::Input) = (;kind=:input,name=value.name,
    recipe=_control_map_recipe_value(value.recipe),states=_control_map_recipe_value(value.states))
_control_map_recipe_value(value::Interaction) = (;kind=:interaction,
    fields=map(_control_map_recipe_value,(; (k=>getfield(value,k) for k in fieldnames(typeof(value)))...)))
function _control_map_recipe_content(recipe::Recipe)
    fields = (; (key=>_control_map_recipe_value(getfield(recipe,key))
        for key in fieldnames(typeof(recipe)) if key!==:scope)...)
    return (;kind=nameof(typeof(recipe)),fields)
end
_control_map_derivation_content(derivation::NamedTuple) = _control_map_recipe_value(derivation)
