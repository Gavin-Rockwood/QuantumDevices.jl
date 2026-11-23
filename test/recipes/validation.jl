@testset "capture owns structural validation for drafts and saved recipes" begin
    d=QuantumDevice(:boundaries)
    q=add_component!(d,make_qubit(:q));r=add_component!(d,make_qubit(:r))
    recipe=Recipe(d,[q])
    @test build(d,recipe).H isa QT.QuantumObject
    for mutate in (
        trial->(trial.states[:q]=[0,0]),
        trial->(trial.states[:missing]=[0]),
        trial->push!(trial.components,:missing),
        trial->push!(trial.inputs.q.owners,:r),
        trial->pop!(trial.inputs.q.operator_names[:q]),
    )
        trial=draft(recipe);mutate(trial)
        @test_throws ArgumentError QuantumDevices.snapshot_recipe(d,trial)
    end
    invalid=draft(recipe);invalid.states[:q]=[0,0]
    @test_throws ArgumentError build(d,invalid)
    @test_throws ArgumentError get_hamiltonian(d,invalid)
    QD=QuantumDevices
    binding=recipe.inputs.q
    badbinding=QD.InputBinding(:other,binding.recipe,binding.owners,recipe.scope,binding.operator_names)
    malformed=QD.CompositeRecipe(recipe.components,(q=badbinding,),recipe.states,recipe.parameters,
        Interaction[],1000,recipe.scope,true,(;))
    @test_throws ArgumentError QD.snapshot_recipe(d,malformed)
    for bad in (
        Combine(inputs=Input[]),
        Combine(inputs=[Input(:same,Leaf(component=:q),[0]),Input(:same,Leaf(component=:r),[0])]),
        Combine(inputs=[Input(:a,Leaf(component=:q),[0]),Input(:b,Leaf(component=:q),[0])]),
        Combine(inputs=[Input(:q,Leaf(component=:q),[0,0])]),
        Combine(inputs=[Input(:q,:absent,[0])]),
        Combine(inputs=[Input(:q,Leaf(component=:q),[0])],interactions=[Interaction(:foreign,op(:r,:x))]),
    )
        @test_throws ArgumentError QD.snapshot_recipe(d,bad)
    end
    d.recipes[:cycle]=Combine(inputs=[Input(:again,:cycle,[0])])
    error=try build(d,:cycle);nothing catch err;err end
    @test error isa ArgumentError
    @test occursin("target",sprint(showerror,error))
    @test occursin("again",sprint(showerror,error))
    saved=add_recipe!(d,:saved,recipe)
    incompatible=QuantumDevice(:other)
    add_component!(incompatible,make_resonator(:q;dimension=3))
    @test_throws ArgumentError build(incompatible,saved)
    @test_throws ArgumentError get_hamiltonian(incompatible,saved)
    @test_throws ArgumentError build(incompatible,GateRecipe(saved;duration=1.))
    @test draft(saved) isa Recipe
    @test ResolvedRecipe(saved,Dict(),Pair{Symbol,ResolvedRecipe}[]).definition isa Recipe
end

@testset "configuration before numerical generator validation and deferred labels" begin
    d=QuantumDevice(:generator_boundary)
    q=add_component!(d,make_qubit(:q))
    calls=Ref(0)
    invalid=Leaf(component=:q,generators=Dict(OperatorRef(:q,:x)=>OperatorGenerator(
        (n;kwargs...)->(calls[]+=1;ones(n,n)),2)))
    @test QuantumDevices.snapshot_recipe(d,invalid) isa Leaf
    @test get_hamiltonian(d,invalid) isa Sym
    @test calls[]==0
    @test_throws ArgumentError build(d,invalid)
    @test calls[]>0
    badname=Leaf(component=:q,generators=Dict(OperatorRef(:q,:unknown)=>OperatorGenerator(n->QT.qeye(n),2)))
    @test_throws ArgumentError QuantumDevices.snapshot_recipe(d,badname)
    @test_throws ArgumentError build(d,Leaf(component=:q,dimension=0))
    labels=Recipe(Recipe(d,[q]);states=[(4,)])
    @test QuantumDevices.snapshot_recipe(d,labels) isa Recipe
    @test get_hamiltonian(d,labels) isa Sym
    @test_throws ArgumentError build(d,labels)
end

@testset "missing parameter use sites remain branch scoped" begin
    d=QuantumDevice(:scoped_missing)
    a=add_component!(d,make_qubit(:a));b=add_component!(d,make_qubit(:b))
    p=Symbolics.variable(:scoped_missing_parameter)
    left=Recipe(d,[a];parameters=Dict(p=>.2))
    right=Recipe(d,[b];interactions=[Interaction(:needs_p,p*b.operators["x"])])
    combined=Recipe((left=left,right=right))
    error=try build(d,combined);nothing catch err;err end
    @test error isa ArgumentError
    @test occursin("needs_p",sprint(showerror,error))
    @test occursin("scoped_missing_parameter",sprint(showerror,error))
    @test build(d,combined;parameters=Dict(p=>.3)).H isa QT.QuantumObject
    rightplain=Recipe(d,[b])
    endpoint=Recipe((left=left,right=rightplain))
    gate=GateRecipe(endpoint;terms=[Interaction(:gate_needs_p,p*b.operators["x"])],duration=1.)
    error=try build(d,gate);nothing catch err;err end
    @test error isa ArgumentError
    @test occursin("gate_needs_p",sprint(showerror,error))
    add_parameters!(d,Dict(p=>0.))
    mixed=GateRecipe(endpoint;terms=[Interaction(:left_p,p*a.operators["x"]),
        Interaction(:right_p,p*b.operators["x"])],duration=1.)
    @test_throws ArgumentError build(d,mixed)
    @test build(d,mixed;parameters=Dict(p=>.3)).parameters[p] == .3
end
