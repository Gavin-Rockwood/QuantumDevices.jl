@testset "recipe input inventory and independent read snapshots" begin
    QD=QuantumDevices
    d=QuantumDevice(:inventory)
    q=add_component!(d,make_qubit(:q;frequency=5.))
    unused=add_component!(d,make_qubit(:unused;frequency=8.))
    ν=q.variables["ν"]
    g=Symbolics.variable(:inventory_g)
    unused_g=Symbolics.variable(:unused_g)
    add_interaction!(d,Interaction(:unused_coupling,unused_g*q.operators["x"]))
    first=Recipe(d,[q];interactions=[Interaction(:shift,g*q.operators["z"])],parameters=Dict(g=>.1))
    second=Recipe(d,[q];parameters=Dict(ν=>5.2))
    add_recipe!(d,:first,first);add_recipe!(d,:second,second)
    @test d.variables[string(g)] === g
    @test d.variables[string(unused_g)] === unused_g
    @test isempty(d.parameters[string(unused.variables["ν"])])
    @test isempty(d.parameters[string(unused_g)])
    @test Set(keys(d.variables))==Set(keys(d.parameters))
    @test Set(keys(d.parameters[string(ν)]))==Set([5.,5.2])
    @test (source=:recipe,path=(:first,),node=()) in d.parameters[string(g)][.1]
    @test (source=:recipe,path=(:first,),node=(:q,)) in d.parameters[string(ν)][5.]
    view=d.parameters
    empty!(view[string(ν)][5.]);empty!(view[string(g)])
    @test !isempty(d.parameters[string(ν)][5.])
    @test haskey(d.parameters[string(g)],.1)
    q.parameters["ν"]=5.1
    @test !haskey(d.parameters[string(ν)],5.1)
    @test haskey(d.parameters[string(ν)],5.)
    @test haskey(d.parameters[string(ν)],5.2)
    missing_recipe=Recipe(d,[q];interactions=[Interaction(:missing,g*q.operators["z"])])
    add_recipe!(d,:missing,missing_recipe)
    @test haskey(d.parameters[string(g)],missing)
    @test_throws ArgumentError build(d,:missing)
    @test build(d,:missing;parameters=Dict(g=>.2)).resolved_parameters[g]==.2
    @test_throws ArgumentError build(d,:first;parameters=Dict(Symbolics.variable(:typo)=>1.))
    @test !isdefined(QuantumDevices,:add_parameters!)
    @test !hasproperty(d,:explicit_parameters)
    @test !hasproperty(Interaction(:plain,g),:parameters)
    @test_throws ArgumentError Interaction(:old,g;parameters=Dict(g=>.1))
end

@testset "inventory is symbolic, declaration complete, and occurrence scoped" begin
    calls=Ref(0)
    p=Symbolics.variable(:declared_input)
    operators,generators=component_operators(:symbolic;x=dimension->(calls[]+=1;error("must not run")))
    component=Component(:symbolic;variables=Dict("p"=>p),parameters=Dict("p"=>2.),
        operators,operator_generators=generators,hamiltonian=val(0)*operators["x"])
    d=QuantumDevice(:symbolic);add_component!(d,component)
    add_recipe!(d,:leaf,Leaf(component=:symbolic))
    @test calls[]==0
    @test haskey(d.parameters[string(p)],2.)
    @test calls[]==0
    @test isnothing(component.dimension)
    child=Recipe(d,:leaf;parameters=Dict(p=>3.))
    parent=Recipe(d,child;parameters=Dict(p=>4.))
    add_recipe!(d,:parent,parent)
    @test length(d.parameters[string(p)][4.])==3
    @test calls[]==0
    # Conflicting numerical values do not hide symbolic declarations.
    shared=Symbolics.variable(:shared_inventory)
    function shared_component(name)
        ops,gens=component_operators(name;z=dimension->QT.sigmaz())
        Component(name;variables=Dict("p"=>shared),parameters=Dict("p"=>1.),operators=ops,
            operator_generators=gens,hamiltonian=shared*ops["z"],dimension=2)
    end
    conflict=QuantumDevice(:conflict)
    a=add_component!(conflict,shared_component(:a));b=add_component!(conflict,shared_component(:b))
    left=Recipe(conflict,[a];parameters=Dict(shared=>2.))
    right=Recipe(conflict,[b];parameters=Dict(shared=>3.))
    conflicting=Recipe(conflict,(left=left,right=right))
    conflict.recipes[:conflicting]=conflicting
    @test conflict.variables[string(shared)] === shared
    error=try conflict.parameters;nothing catch e;e end
    @test error isa ArgumentError
    @test occursin("left",sprint(showerror,error)) && occursin("right",sprint(showerror,error))
end

@testset "gate and historical artifact inventory authority" begin
    QD=QuantumDevices
    d=QuantumDevice(:gates)
    q=add_component!(d,make_qubit(:q;frequency=5.))
    ν=q.variables["ν"]
    base=Recipe(d,[q]);add_recipe!(d,:base,base)
    A=Symbolics.variable(:inventory_drive)
    endpoint=draft(base);endpoint.overrides[A]=0.
    term=Interaction(:drive,A*q.operators["x"])
    gate=GateRecipe(endpoint;terms=[term],controls=[Control(:drive,A,Pulse(Constant(;duration=.01);amplitude=.1))])
    working=QuantumGate(:drive,gate)
    add_gate!(d,(:drive,:small),working)
    @test d.parameters[string(A)][0.] == [(source=:gate,path=(:drive,:small),node=())]
    prepared=build(d,gate;parameters=Dict(A=>.3,ν=>5.3))
    @test prepared.parameters[A]==.3
    reference=QD._gate_reference(prepared.model;parameters=prepared.parameters)
    @test reference.parameters[A]==.3
    problem=CalibrationProblem(prepared,["drive/amplitude"],Matrix{ComplexF64}(I,2,2);
        states=[(0,),(1,)],lb=[0.],ub=[1.])
    result=CalibrationResult(problem,(;u=[.2],retcode=QD.SciMLBase.ReturnCode.MaxIters))
    working.calibration=result
    attach_artifact!(working,:fit,result)
    @test (source=:calibration,path=(:drive,:small),node=()) in d.parameters[string(A)][.3]
    @test (source=:calibration,path=(:drive,:small,:artifacts,:fit),node=()) in d.parameters[string(A)][.3]
    @test (source=:calibration,path=(:drive,:small),node=(:model,)) in d.parameters[string(ν)][5.3]
    add_gate!(d,:idle,QuantumGate(:idle,GateRecipe(base;duration=.01)))
    @test (source=:gate,path=(:idle,),node=()) in d.parameters[string(ν)][5.]
    @test (source=:gate,path=(:idle,),node=(:model,)) in d.parameters[string(ν)][5.]
    attach_artifact!(working,:model,prepared.model)
    @test haskey(d.parameters[string(ν)],5.3)
    @test (source=:model,path=(:drive,:small,:artifacts,:model),node=()) in d.parameters[string(ν)][5.3]
    # A static missing baseline stays missing even when a control replaces it later.
    missing_gate=GateRecipe(base;terms=[term],controls=gate.controls)
    add_gate!(d,:missing,QuantumGate(:missing,missing_gate))
    @test (source=:gate,path=(:missing,),node=()) in d.parameters[string(A)][missing]
    @test_throws ArgumentError build(d,missing_gate)
    @test !haskey(d.recipes[:base].overrides,A)
    @test_throws ArgumentError build(d,endpoint)
end

@testset "derived map inventory uses held physical inputs only" begin
    QD=QuantumDevices
    d=QuantumDevice(:map_inventory)
    q=add_component!(d,make_qubit(:q;frequency=5.))
    ν=q.variables["ν"]
    recipe=Recipe(d,[q])
    context=QD.capture_gate_context(d;include_recipes=false,include_maps=false)
    derivation=(;recipe=QD.snapshot_recipe(d,recipe),context,parameters=QD.RecipeMap(Dict(ν=>5.7)),parameter=ν,
        measurement=()->error("inventory must not measure or invert"))
    map=TabulatedControlMap(:derived;input=:v,output=ν,voltage_knots=(0.,1.),flux_knots=(0.,1.),derivation)
    add_control_map!(d,map)
    usage=d.parameters[string(ν)]
    @test Set(keys(usage))==Set([5.7])
    @test (source=:control_map,path=(:device,:derived),node=()) in usage[5.7]
    @test !haskey(d.variables,"v")
    @test !haskey(usage,0.) && !haskey(usage,1.)
end

@testset "mutated interaction inputs fail symbolic validation" begin
    d=QuantumDevice(:invalid_gate_inputs)
    q=add_component!(d,make_qubit(:q))
    calls=Ref(0)
    for (name,factory) in collect(q.operator_generators)
        q.operator_generators[name]=let f=factory
            dimension -> (calls[]+=1;f(dimension))
        end
    end
    base=Recipe(d,[q])
    A=Symbolics.variable(:invalid_gate_A)
    term=Interaction(:drive,A*q.operators["x"])
    term.variables[string(A)]=Symbolics.variable(:wrong_gate_input)
    gate=GateRecipe(base;terms=[term],duration=.01)
    error=try build(d,gate);nothing catch exception;exception end
    @test error isa ArgumentError
    @test occursin("variables do not match its expression",sprint(showerror,error))
    @test calls[]==0
    add_gate!(d,QuantumGate(:invalid,gate))
    @test_throws ArgumentError d.variables
    @test_throws ArgumentError d.parameters
    @test calls[]==0
end
