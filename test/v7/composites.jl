@testset "explicit composite construction and boundaries" begin
    d = QuantumDevice(:composite)
    q1 = add_component!(d,make_qubit(:q1;frequency=4.7))
    q2 = add_component!(d,make_qubit(:q2;frequency=6.1))
    q3 = add_component!(d,make_qubit(:q3;frequency=7.3))
    pair = Recipe(d,[q1,:q2];states=Dict(q1=>[1,0],:q2=>0:1))
    @test pair isa Recipe
    @test :CompositeRecipe ∉ names(QuantumDevices)
    @test QuantumDevices.CompositeRecipe(d, [q1, q2]) isa Recipe
    @test isempty(d.recipes)
    @test keys(pair.inputs) == (:q1,:q2)
    @test pair.inputs.q1.recipe isa Leaf
    @test propertynames(pair.inputs.q1.operators) |> Set == Set((:x,:y,:z))
    @test :operators in propertynames(pair.inputs.q1)
    @test startswith(sprint(show,pair), "Recipe(")
    g = Symbolics.variable(:local_g)
    add_interaction!(pair,Interaction(:coupling,(g)*pair.inputs.q1.operators.x*pair.inputs.q2.operators.x;parameters=Dict(g=>.13)))
    @test isempty(d.interactions)
    @test_throws ArgumentError add_interaction!(d,Interaction(:bound,pair.inputs.q1.operators.z))
    @test_throws ArgumentError add_interaction!(pair,Interaction(:coupling,pair.inputs.q1.operators.z))
    @test_throws ArgumentError add_interaction!(pair,Interaction(:physical,op(:q1,:z)))
    other = Recipe(d,[q1,q2])
    @test_throws ArgumentError add_interaction!(pair,Interaction(:foreign,other.inputs.q1.operators.z))
    model = build(d,pair)
    reference = 4.7/2*QT.tensor(QT.sigmaz(),QT.qeye(2)) +
        6.1/2*QT.tensor(QT.qeye(2),QT.sigmaz()) + .13*QT.tensor(QT.sigmax(),QT.sigmax())
    @test model.spectrum ≈ real.(QT.eigenenergies(reference))
    @test model.parameters[g] == .13
    @test input_basis_map(model).dims == ([2,2],[4])
    @test all(x -> x isa QT.QuantumObject,values(model.operator_products))
    @test all(x -> x.dims == ([2,2],[2,2]),values(model.operators))
    @test all(r -> r isa OperatorRef,QuantumDevices.refs(model.hamiltonian,:op))
    @test evaluate(model.hamiltonian,model.parameters,Dict(OperatorRef(:q1,:z)=>QT.tensor(QT.sigmaz(),QT.qeye(2)),OperatorRef(:q2,:z)=>QT.tensor(QT.qeye(2),QT.sigmaz()),OperatorRef(:q1,:x)=>QT.tensor(QT.sigmax(),QT.qeye(2)),OperatorRef(:q2,:x)=>QT.tensor(QT.qeye(2),QT.sigmax()))) ≈ reference
    @test build(d,pair;parameters=Dict(g=>0)).spectrum ≈ real.(QT.eigenenergies(reference - .13*QT.tensor(QT.sigmax(),QT.sigmax())))
    @test only(pair.interactions).parameters[string(g)] == .13
    add_recipe!(d,:pair,pair)
    target = Recipe(d,[q1,q2,q3];inputs=(group=:pair,),states=Dict(:group=>[3,0],q3=>[1]))
    @test keys(target.inputs) == (:group,:q3)
    @test propertynames(target.inputs.group.operators) == (:q1,:q2)
    @test target.inputs.group.operators.q2.x isa Sym
    @test_throws ArgumentError target.inputs.group.operators.x
    @test_throws ArgumentError target.inputs.q3.operators.missing
    @test length(build(d,target).spectrum) == 2
    @test build(d,target).recipe.inputs[1].first == :group
    @test length(build(d,Recipe(d,[q1,q2];inputs=(saved=pair,))).spectrum) == 4
    @test length(d.recipes) == 1
    for make in (
        ()->Recipe(d,[q1,q1]),
        ()->Recipe(d,[:absent]),
        ()->Recipe(d,[make_qubit(:q1)]),
        ()->Recipe(d,[q1];inputs=(outside=pair,)),
        ()->Recipe(d,[q1,q2];inputs=(a=pair,b=pair)),
        ()->Recipe(d,[q1,q2,q3];inputs=(q3=pair,)),
        ()->Recipe(d,[q1,q2];inputs=(p=pair,),initial_dims=Dict(q1=>2)),
        ()->Recipe(d,[q1];initial_dims=Dict(q2=>2)),
        ()->Recipe(d,[q1];states=Dict(q2=>[0])),
        ()->Recipe(d,[q1,q2];inputs=(p=pair,),states=Dict(q1=>[0])),
        ()->Recipe(d,[q1];max_dimension=-1),
    )
        @test_throws ArgumentError make()
    end
    for states in (Int[],[0,0],[-1],[0.5])
        @test_throws ArgumentError Recipe(d,[q1];states=Dict(q1=>states))
    end
    @test_throws ArgumentError build(d,Recipe(d,[q1];states=Dict(q1=>[2])))
    @test length(build(d,Recipe(d,[q1];max_dimension=2)).spectrum) == 2
    @test_throws ArgumentError build(d,Recipe(d,[q1];max_dimension=1))
    @test length(build(d,Recipe(d,[q1];max_dimension=0)).spectrum) == 2
    cutoff = add_component!(d,make_resonator(:cutoff;dimension=1000))
    @test QuantumDevices.resolve_node(d,Recipe(d,[cutoff])).dimension == 1000
    @test Recipe(d,[cutoff]).max_dimension == 1000
    @test_throws ArgumentError build(d,Recipe(d,[cutoff];initial_dims=Dict(cutoff=>1001)))
    oversized = Recipe(d, [cutoff]; initial_dims=Dict(cutoff=>1001))
    @test get_hamiltonian(d, oversized) === cutoff.hamiltonian
    oversized_parent = Recipe(d, [cutoff]; inputs=(child=oversized,), states=Dict(:child=>[0]))
    @test get_hamiltonian(d, oversized_parent) === cutoff.hamiltonian
    @test_throws ArgumentError build(d, oversized_parent)
    @test_throws ArgumentError build(d, Combine(inputs=[Input(:q, :pair, [0])], max_dimension=-1))
    huge1 = add_component!(d,make_resonator(:huge1;dimension=typemax(Int)))
    huge2 = add_component!(d,make_resonator(:huge2;dimension=2))
    @test_throws ArgumentError build(d,Recipe(d,[huge1,huge2]))
    error = try build(d,Recipe(d,[huge1,huge2];max_dimension=0));nothing catch e;e end
    @test error isa ArgumentError
    @test occursin(string(big(typemax(Int))*2),sprint(showerror,error))
    # A disabled child limit must not allocate its default retained-state vector
    # before a restrictive ancestor gets the chance to reject it.
    hugechild = Recipe(d,[huge1];max_dimension=0)
    @test_throws ArgumentError build(d,Recipe(d,[huge1];inputs=(huge=hugechild,)))
    d.recipes[:cycle] = Combine(inputs=[Input(:self,:cycle,[0])])
    cycle_error = try build(d,:cycle);nothing catch e;e end
    @test cycle_error isa ArgumentError
    @test occursin("cycle",sprint(showerror,cycle_error))
    @test_throws ArgumentError get_hamiltonian(d,:cycle)
    @test_throws ArgumentError Recipe(d,[q1];inputs=(bad=:cycle,))
    # Named dependencies are checked again at build, including each alias separately.
    add_recipe!(d,:first,Leaf(component=:q1)); add_recipe!(d,:second,Leaf(component=:q2))
    named = Recipe(d,[q1,q2];inputs=(first=:first,second=:second))
    d.recipes[:first] = Leaf(component=:q2)
    d.recipes[:second] = Leaf(component=:q1)
    @test_throws ArgumentError build(d,named)
    # A mutable legacy dependency can introduce a cycle through an object recipe.
    legacy = Combine(inputs=Input[])
    objectcycle = Recipe(d,[q1];inputs=(leaf=Leaf(component=:q1),))
    add_recipe!(d,:objectcycle,objectcycle)
    # Named cycle reaches an object root and must still include a readable path.
    cyclic_legacy = Combine(inputs=[Input(:back,:objectcycle,[0])])
    objectcycle2 = Recipe(d,[q1];inputs=(link=cyclic_legacy,))
    d.recipes[:objectcycle] = objectcycle2
    @test_throws ArgumentError build(d,objectcycle2)
end

@testset "composite parameter inheritance and provenance" begin
    d = QuantumDevice(:parameters)
    q = add_component!(d,make_qubit(:q;frequency=5))
    p = q.variables["ν"]
    leaf = Leaf(component=:q,parameters=Dict(p=>5.1))
    child = Recipe(d,[q];inputs=(q=leaf,),parameters=Dict(p=>5.2))
    localp = Symbolics.variable(:local_only)
    add_interaction!(child,Interaction(:shift,(localp)*child.inputs.q.operators.z;parameters=Dict(localp=>.2)))
    parent = Recipe(d,[q];inputs=(child=child,),parameters=Dict(p=>5.3,localp=>.3))
    model = build(d,parent;parameters=Dict(p=>5.4,localp=>.4))
    @test model.parameters[p] == 5.4
    @test model.parameters[localp] == .4
    @test model.recipe.inputs[1].second.parameters[localp] == .4
    @test !haskey(model.recipe.inputs[1].second.inputs[1].second.parameters,localp)
    @test build(d,parent).parameters[p] == 5.3
    @test build(d,child).parameters[p] == 5.2
    @test build(d,leaf).parameters[p] == 5.1
    @test leaf.parameters[p] == 5.1
    @test q.parameters["ν"] == 5
    @test only(child.interactions).parameters[string(localp)] == .2
    second = Recipe(d,[q];inputs=(child=child,),parameters=Dict(p=>6.3))
    first_model,second_model = build(d,parent),build(d,second)
    @test first_model.parameters[p] == 5.3
    @test second_model.parameters[p] == 6.3
    @test first_model.recipe.inputs[1].second !== second_model.recipe.inputs[1].second
    unrelated = add_component!(d,make_qubit(:unrelated))
    @test !haskey(build(d,parent;parameters=Dict(unrelated.variables["ν"]=>9)).parameters,unrelated.variables["ν"])
    @test_throws ArgumentError build(d,parent;parameters=Dict(Symbolics.variable(:unknown)=>2))
    @test_throws ArgumentError build(d,Recipe(d,[q];parameters=Dict(Symbolics.variable(:unknown)=>2)))
    # Two physical components sharing a parameter must resolve the same value.
    common = Symbolics.variable(:common)
    function shared_component(name)
        Component(name,(frequency=common,),Dict(common=>1.),Dict(OperatorRef(name,:z)=>OperatorGenerator((dim;kwargs...)->QT.sigmaz(),2)),(common)*op(name,:z),"custom",nothing,(;))
    end
    a=add_component!(d,shared_component(:a));b=add_component!(d,shared_component(:b))
    conflicting = Recipe(d,[a,b];inputs=(a=Leaf(component=:a,parameters=Dict(common=>2)),b=Leaf(component=:b,parameters=Dict(common=>3))))
    @test_throws ArgumentError build(d,conflicting)
    @test build(d,conflicting;parameters=Dict(common=>4)).parameters[common] == 4
    # Identical labels in unrelated local scopes are legal.
    ca=Recipe(d,[a]);cb=Recipe(d,[b])
    add_interaction!(ca,Interaction(:same,ca.inputs.a.operators.z))
    add_interaction!(cb,Interaction(:same,cb.inputs.b.operators.z))
    @test length(build(d,Recipe(d,[a,b];inputs=(a=ca,b=cb))).spectrum) == 4
    add_recipe!(d,:qleaf,Leaf(component=:q))
    add_interaction!(d,Interaction(:registered,op(:q,:z)))
    add_recipe!(d,:registered_child,Combine(inputs=[Input(:q,:qleaf,[0,1])],interactions=[:registered]))
    @test_throws ArgumentError build(d,Combine(inputs=[Input(:child,:registered_child,[0,1])],interactions=[:registered]))
    @test_throws ArgumentError build(d,Combine(inputs=[Input(:q,:qleaf,[0,1])],interactions=[:registered,:registered]))
end

@testset "nested complete products preserve staged projection" begin
    d=QuantumDevice(:products)
    r=add_component!(d,make_resonator(:r;dimension=5,frequency=3))
    inner=Recipe(d,[r];states=Dict(r=>[0,1]))
    outer=Recipe(d,[r];inputs=(inner=inner,))
    a=outer.inputs.inner.operators.a
    add_interaction!(outer,Interaction(:product,.3*a*a'))
    model=build(d,outer)
    fulla=QT.destroy(5)
    select=QuantumDevices.selection(5,[1,2])
    complete=select'*(fulla*fulla')*select
    primitive=select'*fulla*select
    expected=QuantumDevices.spectrum_operator([0.,3.])+.3*complete
    incorrect=QuantumDevices.spectrum_operator([0.,3.])+.3*(primitive*primitive')
    @test model.spectrum ≈ real.(QT.eigenenergies(expected))
    @test !(model.spectrum ≈ real.(QT.eigenenergies(incorrect)))
    @test model.operator_products[((OperatorRef(:r,:a),false),)] isa QT.QuantumObject
    # Request a higher ordered word through two intermediate projections.
    middle=Recipe(d,[r];inputs=(inner=inner,),states=Dict(:inner=>[1,0]))
    top=Recipe(d,[r];inputs=(middle=middle,))
    a=top.inputs.middle.operators.a
    add_interaction!(top,Interaction(:quartic,.07*a*a'*a*a'))
    @test build(d,top).spectrum ≈ real.(QT.eigenenergies(QuantumDevices.spectrum_operator([0.,3.])+.07*(select'*(fulla*fulla'*fulla*fulla')*select)))
end

@testset "transmon pair and nested third staged references" begin
    d=QuantumDevice(:three_transmons)
    q1=add_component!(d,make_transmon(:q1,.27,17.;dimension=11))
    q2=add_component!(d,make_transmon(:q2,.29,19.;dimension=11))
    q3=add_component!(d,make_transmon(:q3,.25,16.;dimension=11))
    pair=Recipe(d,[q1,q2];initial_dims=Dict(q1=>11,q2=>11),states=Dict(q1=>0:2,q2=>0:2))
    g12=Symbolics.variable(:g12)
    add_interaction!(pair,Interaction(:coupling12,(g12)*pair.inputs.q1.operators.charge*pair.inputs.q2.operators.charge;parameters=Dict(g12=>.02)))
    leaves=[build(d,Leaf(component=q.name,dimension=11)) for q in (q1,q2,q3)]
    select=QuantumDevices.selection(11,[1,2,3]);eye=QT.qeye(3)
    charges=[select'*input_basis_map(leaf)'*leaf.operators[OperatorRef(q.name,:charge)]*input_basis_map(leaf)*select for (q,leaf) in zip((q1,q2,q3),leaves)]
    energies=[QuantumDevices.spectrum_operator(leaf.spectrum[1:3]) for leaf in leaves]
    expected_pair=QT.tensor(energies[1],eye)+QT.tensor(eye,energies[2])+.02*QT.tensor(charges[1],charges[2])
    pairmodel=build(d,pair)
    @test pairmodel.spectrum ≈ real.(QT.eigenenergies(expected_pair))
    @test input_basis_map(pairmodel)'*expected_pair*input_basis_map(pairmodel) ≈ QuantumDevices.spectrum_operator(pairmodel.spectrum)
    @test pairmodel.operators[OperatorRef(:q2,:charge)] ≈ QT.tensor(eye,charges[2])
    target=Recipe(d,[q1,q2,q3];inputs=(pair=pair,),initial_dims=Dict(q3=>11),states=Dict(:pair=>[0,2,1,4],q3=>0:2))
    g23=Symbolics.variable(:g23)
    add_interaction!(target,Interaction(:coupling23,(g23)*target.inputs.pair.operators.q2.charge*target.inputs.q3.operators.charge;parameters=Dict(g23=>.015)))
    add_recipe!(d,:target,target)
    model=build(d,:target)
    keep=[1,3,2,5];groupselect=QuantumDevices.selection(9,keep)
    groupcharge=groupselect'*(input_basis_map(pairmodel)'*QT.tensor(eye,charges[2])*input_basis_map(pairmodel))*groupselect
    expected=QT.tensor(QuantumDevices.spectrum_operator(pairmodel.spectrum[keep]),eye)+QT.tensor(QT.qeye(4),energies[3])+.015*QT.tensor(groupcharge,charges[3])
    @test model.spectrum ≈ real.(QT.eigenenergies(expected))
    @test input_basis_map(model)'*expected*input_basis_map(model) ≈ QuantumDevices.spectrum_operator(model.spectrum)
    @test input_basis_map(model).dims == ([4,3],[12])
    @test all(x->x isa QT.QuantumObject,values(model.operators))
    @test length(d.recipes) == 1
    @test length(d.interactions) == 0
    @test model.parameters[g12] == .02
    @test model.parameters[g23] == .015
    @test get_hamiltonian(d,:target) isa Sym
end
