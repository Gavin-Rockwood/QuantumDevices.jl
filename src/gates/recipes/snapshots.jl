function capture_gate(device,gate::GateRecipe;frozen=true,components=Dict{Symbol,Component}())
    source = gate.model_recipe isa Symbol ? (isnothing(device) ?
        capture_error(["target"],"supply the device to capture a named gate model") :
        capture_lookup(device.recipes,gate.model_recipe,"recipe",["target"])) : gate.model_recipe
    captured = capture_recipe_tree(device, source, frozen, Any[], ["target"], components)
    model = captured.recipe
    owners = captured.owners
    operators = captured.operators
    old_scope = source isa CompositeRecipe ? source.scope : nothing
    new_scope = model isa CompositeRecipe ? model.scope : nothing
    extra = [copy_interaction(term,frozen;old_scope,new_scope) for term in gate.terms]
    for term in extra
        capture_operator_refs(term,model,owners,operators,["target","gate $(term.name)"];local_scope=false)
    end
    controls = deepcopy(collect(gate.controls))
    return _physical_gate(model,frozen ? Tuple(controls) : controls,
        frozen ? RecipeVector(extra) : extra,gate.duration,
        gate.duration_policy,
        frozen ? freeze_recipe_value(gate.metadata) : thaw_recipe_value(gate.metadata))
end

draft(gate::GateRecipe) = capture_gate(nothing, gate; frozen=false)


# Physical definitions are provenance, separate from mutable numerical projections.
# Callable factories retain their bindings; external mutable closure state cannot
# be made reproducible by copying a device.
function capture_gate_context(device; include_recipes=true, include_maps=true)
    components = Dict(name => (; name=c.name, variables=RecipeMap(c.variables), parameters=RecipeMap(c.parameters),
        operators=RecipeMap(c.operators),
        operator_generators=RecipeMap(Dict(k=>OperatorGenerator(deepcopy(v.f),v.default_dimension)
            for (k,v) in c.operator_generators)), hamiltonian=deepcopy(c.hamiltonian),
        type=c.type, label=c.label, metadata=freeze_recipe_value(c.metadata))
        for (name,c) in device.components)
    return (; name=device.name, components=RecipeMap(components),
        explicit_parameters=RecipeMap(device.explicit_parameters),
        control_maps=RecipeMap(include_maps ? device.control_maps : Dict{Symbol,AbstractControlMap}()),
        interactions=RecipeMap(Dict(name=>copy_interaction(i,true) for (name,i) in device.interactions)),
        recipes=RecipeMap(Dict(name=>(is_recipe_snapshot(r) ? snapshot_recipe(nothing,r) : copy_recipe_tree(device,r,true))
            for (name,r) in (include_recipes ? device.recipes : Dict{Symbol,Recipe}()))),
        metadata=freeze_recipe_value(device.metadata))
end
