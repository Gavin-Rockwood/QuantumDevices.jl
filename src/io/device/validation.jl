# Validate stored coordinates and dimensions without generating operators or propagating states.
function _validate_stored_model(model,path)
    n=size(model.H,1)
    size(model.H)==(n,n) && length(model.spectrum)==n && length(model.states)==n || _ioerror(path,"inconsistent model dimensions")
    for state in model.states
        length(state.data)==n && state.dimensions.to==model.H.dimensions.to || _ioerror(path,"state basis differs from stored H")
    end
    for operator in Iterators.flatten((values(model.operators),values(model.operator_products)))
        size(operator)==(n,n) && operator.dimensions==model.H.dimensions || _ioerror(path,"cached operator basis differs from stored H")
    end
    validate_parameter_values(model.parameters)
    return model
end
function _validate_stored_prepared(prepared,path)
    _validate_stored_model(prepared.model,path*"/model")
    validate_parameter_values(prepared.parameters)
    for control in prepared.recipe.controls
        haskey(prepared.parameters,control.target) || _ioerror(path,"control target has no captured physical value")
        for (_,(_,value)) in calibration_parameters(control.pulse)
            value isa Num && !issubset(Set(Symbolics.get_variables(value)),Set(keys(prepared.parameters))) && _ioerror(path,"offset references absent captured values")
        end
    end
    !isnothing(prepared.context) && _validate_device_definitions(QuantumDevice(prepared.context)) # Declaration/index validation only.
    return prepared
end
function _validate_stored_evolution(evolution,path)
    n=length(evolution.reference.energies)
    m=isnothing(evolution.inputs) ? n : length(evolution.inputs)
    length(evolution.times)==length(evolution.propagators)>0 || _ioerror(path,"propagation samples do not align")
    all(t->t isa Real && isfinite(t) && 0<=t<=evolution.recipe.duration,evolution.times) &&
        issorted(evolution.times) && allunique(evolution.times) && last(evolution.times)==evolution.recipe.duration || _ioerror(path,"invalid saved evolution times")
    all(U->size(U)==(n,m),evolution.propagators) || _ioerror(path,"invalid propagation dimensions")
    return evolution
end
function _validate_stored_result(result,path)
    _validate_stored_prepared(result.provenance.prepared,path*"/prepared")
    count=length(result.coordinates)
    count==length(result.fitted_values)==length(result.solution.u) && allunique(result.coordinates) || _ioerror(path,"invalid optimizer coordinate count")
    all(_gate_real,result.fitted_values) || _ioerror(path,"fitted coordinates must be finite real values")
    settings=result.settings.problem
    for bound in (settings.lb,settings.ub)
        isnothing(bound) || length(bound)==count || _ioerror(path,"bounds do not match coordinates")
    end
    for i in 1:count
        lower=isnothing(settings.lb) ? -Inf : settings.lb[i]
        upper=isnothing(settings.ub) ? Inf : settings.ub[i]
        lower<=result.fitted_values[i]<=upper || _ioerror(path,"stored fitted values outside bounds")
    end
    length(result.states)==size(result.target,1)==size(result.target,2) || _ioerror(path,"target and ordered states differ")
    !isnothing(result.state_keys) && length(result.state_keys)!=length(result.states) && _ioerror(path,"state labels differ from ordered kets")
    n=size(result.provenance.prepared.model.H,1)
    all(s->length(s.data)==n,result.states) || _ioerror(path,"stored input ket dimension differs")
    if result.parameterization.mode===:symbolic
        all(is_parameter,result.coordinates) || _ioerror(path,"nonvariable symbolic coordinates")
        _check_calibration_symbols(result.provenance.prepared,result.coordinates)
        for (target,expression) in result.parameterization.links
            _gate_numeric_leaf(result.recipe.controls,target)
            issubset(Set(Symbolics.get_variables(expression)),Set(result.coordinates)) || _ioerror(path,"link references undeclared optimization variables")
        end
    elseif result.parameterization.mode===:direct
        Tuple(result.coordinates)==result.paths && isempty(result.parameterization.links) || _ioerror(path,"invalid direct parameterization")
    else
        _ioerror(path,"unknown parameterization mode")
    end
    if !isnothing(result.validation) && !isnothing(result.validation.evolution)
        _validate_stored_evolution(result.validation.evolution,path*"/validation/evolution")
    end
    return result
end
function _validate_stored_gate(gate,path)
    !isnothing(gate.calibration) && _validate_stored_result(gate.calibration,path*"/calibration")
    for (name,artifact) in gate.artifacts
        artifact isa NumericalModel && _validate_stored_model(artifact,"$path/artifacts/$name")
        artifact isa CalibrationResult && _validate_stored_result(artifact,"$path/artifacts/$name")
        artifact isa GateEvolution && _validate_stored_evolution(artifact,"$path/artifacts/$name")
        artifact isa Union{NumericalModel,CalibrationResult,GateEvolution,TrackingResult,DressedStates,_AnalysisArtifact} || _ioerror(path,"unsupported attached artifact")
    end
    return gate
end

function _validate_recipe_context(recipe,device,path,active=Base.IdSet())
    if recipe isa Symbol
        haskey(device.recipes,recipe) || _ioerror(path,"unknown named recipe $recipe")
        return _validate_recipe_context(device.recipes[recipe],device,path*"/$recipe",active)
    end
    recipe in active && _ioerror(path,"cyclic recipe definitions")
    push!(active,recipe)
    try
        if recipe isa Leaf
            haskey(device.components,recipe.component) || _ioerror(path,"unknown physical component $(recipe.component)")
        elseif recipe isa Combine
            allunique(input.name for input in recipe.inputs) || _ioerror(path,"duplicate input names")
            for input in recipe.inputs;_validate_recipe_context(input.recipe,device,path*"/$(input.name)",active);end
        elseif recipe isa CompositeRecipe
            all(name->haskey(device.components,name),recipe.components) || _ioerror(path,"unknown physical owner")
            for input in values(recipe.inputs);_validate_recipe_context(input.recipe,device,path*"/$(input.name)",active);end
        else
            _ioerror(path,"unsupported recipe representation")
        end
    finally
        delete!(active,recipe)
    end
end
function _validate_device_definitions(device)
    validate_default_sources(device_default_sources(device))
    for (name,recipe) in device.recipes
        _validate_recipe_context(recipe,device,"recipes/$name")
    end
    function visit(repository,path)
        for (name,node) in repository
            node isa Dict ? visit(node,"$path/$name") : _validate_recipe_context(node.recipe.model_recipe,device,"$path/$name/model_recipe")
        end
    end
    _validate_gate_repository(device.gates)
    visit(device.gates,"gates")
    return device
end
