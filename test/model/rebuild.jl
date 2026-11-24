struct CountedProducts
    data::Dict{Tuple,QT.QuantumObject}
    reads::Dict{Tuple,Int}
end
function Base.getindex(products::CountedProducts,word)
    products.reads[word] = get(products.reads,word,0)+1
    return products.data[word]
end

@testset "stable model preparation, stored H and shared gate operating points" begin
    d = QuantumDevice(:rebuild)
    q = add_component!(d,make_qubit(:q;frequency=4.))
    base = Recipe(d,[q])
    reduction = Recipe(base;states=[(1,),(0,)])
    captured = QuantumDevices.snapshot_recipe(d,reduction)
    node,_ = QuantumDevices.prepare_node(d,captured,Dict())
    requests = deepcopy(node.selections)
    required = copy(node.required_products)
    physicalH = node.hamiltonian
    firstmodel = QuantumDevices.assemble_node(node)
    secondmodel = QuantumDevices.assemble_node(node)
    @test isequal(node.selections,requests)
    @test node.selections == Any[[(1,),(0,)]]
    @test node.required_products == required
    @test node.hamiltonian === physicalH
    @test firstmodel.H === firstmodel.H
    @test make_hamiltonian(firstmodel) !== firstmodel.H
    @test make_hamiltonian(firstmodel).data == firstmodel.H.data
    @test secondmodel.H ≈ firstmodel.H
    # Pure reduction needs no H/ramp inputs: these sentinels would make any
    # attempted solve/interaction preparation fail. Upstream eigenpairs suffice.
    inherited = QuantumDevices.composite_eigensystem(reduction,[build(d,base)],[[1,2]],
        nothing,nothing,nothing)
    @test inherited.energies[:] ≈ firstmodel.energies[:]
    @test inherited.tracking.total_steps == build(d,base).labeling.total_steps
    @test count(==(:H),propertynames(firstmodel)) == 1
    @test !hasproperty(node,:requests)
    @test_throws MethodError NumericalModel(firstmodel.recipe,collect(firstmodel.spectrum),
        firstmodel.operators,firstmodel.hamiltonian,firstmodel.owners,input_basis_map(firstmodel),firstmodel.operator_products)
    converted = make_hamiltonian(firstmodel)
    converted.data[1,1] += .7
    @test converted.data != firstmodel.H.data
    # Gate-only declarations/operating points share precedence with ordinary H.
    p = Symbolics.variable(:rebuild_extra)
    upstream = Recipe(d,[q];parameters=Dict(p=>.2))
    downstream = Recipe(upstream;states=0:1,parameters=Dict(p=>.3))
    term = Interaction(:extra,p*q.operators["x"])
    gate = GateRecipe(downstream;terms=[term],duration=1.)
    prepared = build(d,gate)
    @test prepared.parameters[p] == .3
    @test build(d,gate;parameters=Dict(p=>.4)).parameters[p] == .4
    @test prepared.model.H ≈ build(d,Recipe(base;states=0:1)).H
    @test make_hamiltonian(prepared) ≈ prepared.model.H+.3*prepared.model.operators[OperatorRef(:q,:x)]
    @test make_hamiltonian(prepared.model,gate) ≈ make_hamiltonian(prepared)
    upstreamgate = GateRecipe(Recipe(upstream;states=0:1);terms=[term],duration=1.)
    @test build(d,upstreamgate).parameters[p] == .2
    @test make_hamiltonian(build(d,upstreamgate).model,upstreamgate) ≈ make_hamiltonian(build(d,upstreamgate))
    changed = draft(downstream)
    changed.overrides[p] = 1.
    changed_gate = GateRecipe(changed;terms=[Interaction(:changed,p*q.operators["x"])],duration=1.)
    @test build(d,changed_gate).parameters[p] == 1.
    @test build(d,changed_gate;parameters=Dict(p=>.4)).parameters[p] == .4
end

@testset "complete local products projected once and stored reference cleanup" begin
    ax,ay,bx,by = OperatorRef(:a,:x),OperatorRef(:a,:y),OperatorRef(:b,:x),OperatorRef(:b,:y)
    A = QT.Qobj([0. 1 0;1 0 2;0 2 1])
    B = QT.Qobj([1. 0 1;0 2 0;1 0 0])
    leftword = ((ax,false),(ay,false))
    rightx,righty = ((bx,false),),((by,false),)
    left = CountedProducts(Dict{Tuple,QT.QuantumObject}(leftword=>A*B),Dict{Tuple,Int}())
    right = CountedProducts(Dict{Tuple,QT.QuantumObject}(rightx=>QT.sigmax(),righty=>QT.sigmay()),Dict{Tuple,Int}())
    required = Set{Tuple}([(leftword...,rightx...),(leftword...,righty...),leftword,rightx])
    P = QuantumDevices.selection(3,[1,2])
    models = [(;owners=Set([:a]),operator_products=left),(;owners=Set([:b]),operator_products=right)]
    products = QuantumDevices.projected_products(models,[P,QT.qeye(2)],required)
    @test all(==(1),values(left.reads))
    @test all(==(1),values(right.reads))
    @test Set(keys(left.reads)) == Set([leftword])
    @test Set(keys(right.reads)) == Set([rightx,righty])
    @test products[(leftword...,rightx...)] ≈ QT.tensor(P'*A*B*P,QT.sigmax())
    @test norm(products[(leftword...,rightx...)]-QT.tensor((P'*A*P)*(P'*B*P),QT.sigmax()))>1e-4
    d = QuantumDevice(:stored_cleanup)
    raw = QT.Qobj(ComplexF64[1+1e-11im 0;0 2])
    operators,generators = component_operators(:near;x=(n->raw))
    component = add_component!(d,Component(:near;variables=Dict(),parameters=Dict(),operators,
        operator_generators=generators, dimension=2,hamiltonian=operators["x"]))
    model = build(d,Leaf(component=:near))
    @test ishermitian(model.H)
    @test make_hamiltonian(model) == model.H
    @test make_hamiltonian(model;herm_check=0).data == raw.data
    @test model.operators[OperatorRef(:near,:x)].data == raw.data
    @test all(norm(model.H*ket-energy*ket)<1e-12 for (ket,energy) in zip(model.states,model.energies))
end
