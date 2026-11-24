# Frozen wrappers are ownership, not part of a physical definition fingerprint.
function _canonical_recipe_containers(value)
    if value isa AbstractDict
        result=Dict(k=>_canonical_recipe_containers(v) for (k,v) in value)
        get(result,"tag",nothing)=="recipe_map" && (result["tag"]="dict")
        return result
    end
    value isa AbstractVector && return [_canonical_recipe_containers(v) for v in value]
    return value
end
function canonical_recipe_record(recipe,w,path)
    fields=Dict{String,Any}()
    for field in fieldnames(typeof(recipe))
        value=getfield(recipe,field)
        fields[String(field)]=field===:scope ? Dict("tag"=>"scope","id"=>_scope(w,value)) : _encode_record(value,w,"$path/$field")
    end
    return _canonical_recipe_containers(Dict("kind"=>string(nameof(typeof(recipe))),"fields"=>fields))
end
function intern_recipe!(w,recipe,path)
    recipe in w.active && _ioerror(path,"cyclic recipe references")
    push!(w.active,recipe)
    parent=w.scopes
    w.scopes=Dict{Symbol,String}()
    local_scopes=nothing
    record=nothing
    try
        record=canonical_recipe_record(recipe,w,path)
        local_scopes=w.scopes
    finally
        w.scopes=parent
        delete!(w.active,recipe)
    end
    id=recipe_fingerprint(record)
    haskey(w.recipes,id) && _canonical(w.recipes[id])!=_canonical(record) && _ioerror(path,"recipe digest collision")
    w.recipes[id]=record
    aliases=Dict(local_id=>_scope(w,scope) for (scope,local_id) in sort!(collect(local_scopes);by=last))
    return Dict("tag"=>"recipe_ref","id"=>id,"scopes"=>aliases)
end
function _restore_recipe_reference(r,reference,path)
    id=reference["id"]
    haskey(r.recipes,id) || _ioerror(path,"unknown recipe node $id")
    id in r.active && _ioerror(path,"cyclic recipe nodes")
    push!(r.active,id)
    parent=r.scopes
    r.scopes=Dict(local_id=>_unScope(r,target) for (local_id,target) in reference["scopes"])
    try
        record=r.recipes[id]
        fields=record["fields"]
        types=Dict("Leaf"=>Leaf,"Combine"=>Combine,"CompositeRecipe"=>CompositeRecipe)
        haskey(types,record["kind"]) || _ioerror(path,"unknown recipe kind")
        Set(keys(fields))==Set(String.(fieldnames(types[record["kind"]]))) || (record["kind"]=="Leaf" && Set(keys(fields))==setdiff(Set(String.(fieldnames(Leaf))),Set(["input_defaults"]))) || _ioerror(path,"recipe field schema mismatch")
        values=Dict(Symbol(k)=>(k=="scope" ? _unScope(r,v["id"]) : _decode_record(v,r,"$path/$k")) for (k,v) in fields)
        kind=record["kind"]
        kind=="Leaf" && return Leaf(;component=values[:component],dimension=values[:dimension],generators=_io_owned_map(values[:generators]),parameters=_io_owned_map(values[:parameters]),control_maps=_io_owned_map(values[:control_maps]),input_defaults=isnothing(get(values,:input_defaults,nothing)) ? nothing : _io_owned_map(values[:input_defaults]))
        kind=="Combine" && return Combine(_io_owned_vector(values[:inputs]),_io_owned_vector(values[:interactions]),_io_owned_map(values[:parameters]),values[:max_dimension],values[:dressing],values[:dressing_options],_io_owned_map(values[:control_maps]))
        kind=="CompositeRecipe" && return CompositeRecipe(_io_owned_vector(values[:components]),values[:inputs],_io_owned_map(values[:states]),_io_owned_map(values[:parameters]),_io_owned_vector(values[:interactions]),values[:max_dimension],values[:scope],values[:dressing],values[:dressing_options],_io_owned_map(values[:control_maps]))
        _ioerror(path,"unknown recipe kind $kind")
    finally
        r.scopes=parent
        delete!(r.active,id)
    end
end

# Working gates recover editable containers from this independently decoded owner.
_io_working(x::RecipeMap)=Dict(_io_working(k)=>_io_working(v) for (k,v) in x)
_io_working(x::RecipeVector)=[_io_working(v) for v in x]
_io_working(x::NamedTuple)=map(_io_working,x)
_io_working(x::Tuple)=map(_io_working,x)
_io_working(x)=x
function _io_working(x::Leaf)
    Leaf(component=x.component,dimension=x.dimension,generators=_io_working(x.generators),parameters=_io_working(x.overrides),control_maps=_io_working(x.control_maps),input_defaults=x.input_defaults)
end
function _io_working(x::InputBinding)
    InputBinding(x.name,_io_working(x.recipe),_io_working(x.owners),x.scope,_io_working(x.operator_names))
end
_io_working(x::Input)=Input(x.name,_io_working(x.recipe),_io_working(x.states))
function _io_working(x::CompositeRecipe)
    CompositeRecipe(_io_working(x.components),_io_working(x.inputs),_io_working(x.states),_io_working(x.overrides),_io_working(x.interactions),x.max_dimension,x.scope,x.dressing,_io_working(x.dressing_options),_io_working(x.control_maps))
end
_io_working(x::Combine)=Combine(_io_working(x.inputs),_io_working(x.interactions),_io_working(x.overrides),x.max_dimension,x.dressing,_io_working(x.dressing_options),_io_working(x.control_maps))
function _io_working(x::GateRecipe)
    _physical_gate(_io_working(x.model_recipe),collect(x.controls),collect(Interaction,x.terms),x.duration,x.duration_policy,_io_working(x.metadata))
end

function restore_recipe_reference(r,reference,path;ownership=:snapshot)
    ownership in (:snapshot,:working) || _ioerror(path,"unknown recipe ownership mode")
    recipe=_restore_recipe_reference(r,reference,path)
    return ownership===:working ? _io_working(recipe) : recipe
end
