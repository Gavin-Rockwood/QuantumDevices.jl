# Both modes expose native coordinates and actual written pulse paths separately.
function _direct_parameterization(gate::GateRecipe,paths)
    paths isa AbstractString && throw(ArgumentError("pass an ordered collection of pulse paths"))
    names = Tuple(String.(collect(paths)))
    !isempty(names) || throw(ArgumentError("select nonempty unique pulse paths"))
    normalized = map(names) do path
        _,canonical,_,value = _gate_numeric_leaf(gate.controls,path)
        return canonical=>value
    end
    coordinates = Tuple(first.(normalized))
    allunique(coordinates) || throw(ArgumentError("duplicate pulse paths or aliases"))
    return (;coordinates,paths=coordinates,
        parameterization=(;mode=:direct,links=RecipeMap(Dict{String,Any}())),
        initial_values=Tuple(last.(normalized)))
end

# Symbolic optimization coordinates never resolve against physical defaults.
function _check_calibration_symbols(prepared::PreparedGate,coordinates::Tuple)
    physical = Set(keys(prepared.parameters))
    for expression in vcat([prepared.model.hamiltonian],prepared.terms)
        union!(physical,refs(expression,:param))
    end
    for control in prepared.recipe.controls, (_,(_,value)) in parameters(control.pulse)
        value isa Num && push!(physical,value)
    end
    isempty(intersect(Set(coordinates),physical)) || throw(ArgumentError(
        "optimization variables must be distinct from physical model, term and offset references"))
    return nothing
end

function _symbolic_parameterization(prepared::PreparedGate,variables,links)
    variables isa Union{AbstractVector,Tuple} || throw(ArgumentError(
        "variables must be an explicitly ordered vector or tuple of variable=>initial pairs; convert a Dict to an ordered pair list"))
    !isempty(variables) && all(p -> p isa Pair,variables) || throw(ArgumentError("declare nonempty variable=>initial pairs"))
    coordinates = Tuple(first.(variables))
    all(is_parameter,coordinates) && allunique(coordinates) || throw(ArgumentError("optimization variables must be unique atomic real Symbolics variables"))
    initial_values = Tuple(last.(variables))
    all(_gate_real,initial_values) || throw(ArgumentError("initial coordinates must be finite concrete real values"))
    links isa AbstractDict && !isempty(links) || throw(ArgumentError("links must be a nonempty dictionary of physical pulse paths to scalar expressions"))
    normalized = Dict{String,Any}()
    for (path,expression) in links
        path isa String || throw(ArgumentError("link paths must be Strings"))
        _,canonical,_,_ = _gate_numeric_leaf(prepared.recipe.controls,path)
        haskey(normalized,canonical) && throw(ArgumentError("duplicate link target through aliases: $path"))
        expression isa Real && !(expression isa Bool) || throw(ArgumentError("link $path must be a real scalar expression"))
        normalized[canonical] = expression
    end
    _check_calibration_symbols(prepared,coordinates)
    parameterization = (;mode=:symbolic,links=RecipeMap(normalized))
    _parameter_updates(coordinates,parameterization,initial_values)
    return (;coordinates,paths=Tuple(sort!(collect(keys(normalized)))),parameterization,initial_values)
end

function _parameter_updates(coordinates::Tuple,parameterization::NamedTuple,values)
    length(coordinates)==length(values) || throw(DimensionMismatch("wrong coordinate count"))
    all(_gate_real,values) || throw(ArgumentError("optimizer coordinates must be finite concrete real values"))
    if parameterization.mode === :direct
        return [path=>value for (path,value) in zip(coordinates,values)]
    end
    parameterization.mode === :symbolic || throw(ArgumentError("unknown calibration parameterization mode"))
    substitutions = Dict(variable=>value for (variable,value) in zip(coordinates,values))
    updates = Pair{String,Any}[]
    for (path,expression) in parameterization.links
        value = try
            scalar_value(expression,substitutions)
        catch error
            error isa InterruptException && rethrow()
            throw(ArgumentError("cannot resolve link $path: $(sprint(showerror,error))"))
        end
        _gate_real(value) || throw(ArgumentError("link $path must resolve to a finite concrete real value"))
        push!(updates,path=>value)
    end
    return updates
end
