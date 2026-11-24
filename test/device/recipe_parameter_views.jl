@testset "captured recipe and model input views" begin
    d=QuantumDevice(:captured_views)
    q=add_component!(d,make_qubit(:q;frequency=5.))
    ν=q.variables["ν"]
    recipe=Recipe(d,[q])
    @test recipe.variables[string(ν)] === ν
    @test Set(keys(recipe.variables))==Set(keys(recipe.parameters))
    @test Set(keys(recipe.parameters[string(ν)]))==Set([5.])
    @test (source=:recipe,path=(),node=(:q,)) in recipe.parameters[string(ν)][5.]
    @test (source=:recipe,path=(),node=()) in recipe.parameters[string(ν)][5.]
    q.parameters["ν"]=6.
    @test haskey(recipe.parameters[string(ν)],5.)
    model=build(d,recipe)
    @test model.resolved_parameters[ν]==5.
    @test model.recipe.resolved_parameters[ν]==5.
    @test isequal(model.variables,model.recipe.variables)
    @test isequal(model.variables,recipe.variables)
    @test model.parameters==model.recipe.parameters
    @test haskey(model.parameters[string(ν)],5.)
    symbolic_hamiltonian=make_hamiltonian(model,[ν];auto_time_dependence=false)
    named_hamiltonian=make_hamiltonian(model,[string(ν)];auto_time_dependence=false)
    runtime=Dict(ν=>5.25)
    @test named_hamiltonian(Dict(string(ν)=>5.25),0.).data ≈ symbolic_hamiltonian(runtime,0.).data
    @test_throws KeyError named_hamiltonian(runtime,0.)
    @test_throws ArgumentError make_hamiltonian(model,[ν,string(ν)])
    @test_throws ArgumentError make_hamiltonian(model,[string(ν),ν])
    @test_throws ArgumentError make_hamiltonian(model,[ν,ν])
    dynamic_named=make_hamiltonian(model,[string(ν)])
    dynamic_symbolic=make_hamiltonian(model,[ν])
    @test dynamic_named(Dict(string(ν)=>(t->5.25+t)),.2) ≈ dynamic_symbolic(Dict(ν=>(t->5.25+t)),.2)
    @test_throws ArgumentError make_hamiltonian(model,["unknown_input"])
    @test_throws ArgumentError make_hamiltonian(model,string(ν))
    @test haskey(build(d,recipe;parameters=Dict(ν=>7.)).parameters[string(ν)],7.)
    recipe.overrides[ν]=5.5
    @test haskey(recipe.parameters[string(ν)],5.5)
    @test model.resolved_parameters[ν]==5.
    copy=recipe.parameters
    empty!(copy[string(ν)][5.5]);empty!(copy)
    @test !isempty(recipe.parameters[string(ν)][5.5])
    saved=add_recipe!(d,:saved,recipe)
    @test_throws MethodError setindex!(d.recipes[:saved].overrides,8.,ν)
    editable=draft(d.recipes[:saved])
    editable.overrides[ν]=8.
    @test haskey(editable.parameters[string(ν)],8.)
    @test haskey(d.recipes[:saved].parameters[string(ν)],5.5)
    unbound=Leaf(component=:q)
    @test_throws ArgumentError unbound.variables
    @test_throws ArgumentError unbound.parameters
    @test build(d,unbound).resolved_parameters[ν]==6.
    # A known empty declaration set is distinct from an unbound leaf.
    empty_leaf=Leaf(component=:empty,input_defaults=QuantumDevices.RecipeMap(Dict{Symbolics.Num,Number}()))
    @test isempty(empty_leaf.variables)
    @test isempty(empty_leaf.parameters)
    g=Symbolics.variable(:mixed_runtime_g)
    coupled=build(d,Recipe(d,[q];interactions=[Interaction(:shift,g*q.operators["z"])],parameters=Dict(g=>.1)))
    mixed_static=make_hamiltonian(coupled,[string(ν),g];auto_time_dependence=false)
    symbols_static=make_hamiltonian(coupled,[ν,g];auto_time_dependence=false)
    @test mixed_static(Dict{Any,Any}(string(ν)=>6.2,g=>.3),.2) ≈ symbols_static(Dict(ν=>6.2,g=>.3),.2)
    mixed_dynamic=make_hamiltonian(coupled,[string(ν),g])
    symbols_dynamic=make_hamiltonian(coupled,[ν,g])
    @test mixed_dynamic(Dict{Any,Any}(string(ν)=>(t->6.0 + t),g=>(t->.1+t)),.2) ≈
        symbols_dynamic(Dict(ν=>(t->6.0 + t),g=>(t->.1+t)),.2)
    @test_throws ArgumentError make_hamiltonian(coupled,[g,string(g)])
    q.variables["ν"]=Symbolics.variable(:changed_input)
    @test_throws ArgumentError build(d,recipe)
end
