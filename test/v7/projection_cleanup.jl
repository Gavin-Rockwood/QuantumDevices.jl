@testset "ordered selections and deferred tracking lookup" begin
    d=QuantumDevice(:projection_selections)
    a=add_component!(d,make_resonator(:a;frequency=2.,dimension=4))
    b=add_component!(d,make_resonator(:b;frequency=3.,dimension=4))
    left=Recipe(d,[a]);right=Recipe(d,[b])
    lm,rm=build(d,left),build(d,right)
    @test QuantumDevices.resolve_selection(lm,[3,1],:left)==[4,2]
    @test QuantumDevices.resolve_selection(lm,[(3,),(1,)],:left)==[4,2]
    @test QuantumDevices.resolve_selection(lm,nothing,:left)==[1,2,3,4]
    for (request,fragment) in (([(9,)],"unavailable"),([(0,0)],"requires 1 levels"))
        error=try QuantumDevices.resolve_selection(lm,request,:left);nothing catch err;err end
        @test error isa ArgumentError
        @test occursin("left",sprint(showerror,error))
        @test occursin(fragment,sprint(showerror,error))
    end
    untracked=build(d,Recipe(d,[a];dressing=false))
    @test_throws ArgumentError QuantumDevices.resolve_selection(untracked,[(0,)],:untracked)
    @test QuantumDevices.resolve_selection(untracked,[2,0],:untracked)==[3,1]
    recipe=Recipe((left=left,right=right);states=Dict(:left=>[3,1],:right=>[2,0]))
    node,_=QuantumDevices.prepare_node(d,recipe,Dict())
    requests=deepcopy(node.selections)
    indices,projectors=QuantumDevices.selected_inputs(node,[lm,rm])
    @test indices==[[4,2],[3,1]]
    @test isequal(node.selections,requests)
    @test projectors[1].dims==([4],[2])
    @test projectors[1]*QT.basis(2,0) ≈ lm.states[4]
    @test projectors[1]*QT.basis(2,1) ≈ lm.states[2]
    model=build(d,recipe)
    @test Matrix(model.H.data) ≈ Diagonal([12.,6.,8.,2.])
    @test Matrix(model.operators[OperatorRef(:a,:n)].data) ≈ Diagonal([3.,3.,1.,1.])
    @test model.H.dims==([2,2],[2,2])
    for (energy,ket) in zip(model.spectrum,model.states)
        @test model.H*ket ≈ energy*ket
    end
end

@testset "stable degenerate reduction order and inherited diagnostics" begin
    d=QuantumDevice(:degenerate_projection)
    a=add_component!(d,make_qubit(:a;frequency=1.))
    b=add_component!(d,make_qubit(:b;frequency=1.))
    base=Recipe(d,[a,b])
    model=build(d,base)
    source=model.eigensystem.tracking
    minimums=[.99,.8,.7,.6]
    logs=log.(minimums)
    accuracy=merge(getfield(source,:_accuracy),
        (minimum_overlaps=minimums,cumulative_log_overlaps=logs))
    diagnostics=[LabelDiagnostic(minimums[i],logs[i],false,source.diagnostics[i].degenerate) for i in 1:4]
    tracking=TrackingResult(source.keys,source.saved_steps,source.total_steps,collect(model.states),
        (energies=collect(model.spectrum),),collect(0:3),accuracy,diagnostics)
    decorated=NumericalModel(model.recipe,model.H,DressedStates(tracking,model.eigensystem.components),
        model.operators,model.hamiltonian,model.owners,model.operator_products)
    reduction=Recipe(base;states=[2,1,3],dressing_options=(trajectory=x->error("no new tracking"),))
    node,_=QuantumDevices.prepare_node(d,reduction,Dict())
    indices,projectors=QuantumDevices.selected_inputs(node,[decorated])
    @test indices==[[3,2,4]] # equal-energy rows keep request order
    @test projectors[1].dims==([2,2],[3])
    @test projectors[1].data==hcat((decorated.states[i].data for i in [3,2,4])...)
    # No H/ramp inputs exist here: a solve or tracking pass would fail.
    reduced=QuantumDevices.composite_eigensystem(node.definition,[decorated],indices,nothing,nothing,nothing)
    @test reduced.keys==Tuple(source.keys[i] for i in [3,2,4])
    @test reduced.energies[:] == model.spectrum[[3,2,4]]
    @test reduced.tracking.diagnostics[:] == diagnostics[[3,2,4]]
    @test reduced.tracking.accuracy.minimum_overlaps[:] == minimums[[3,2,4]]
    @test reduced.tracking.accuracy.cumulative_log_overlaps[:] == logs[[3,2,4]]
    @test all(isapprox.(collect(reduced.states),[QT.basis(3,i) for i in 0:2]))
    @test reduced.tracking.total_steps==source.total_steps
    @test source.keys==model.eigensystem.keys
    @test getfield(source,:_accuracy).minimum_overlaps==ones(4)
    actual=build(d,reduction)
    @test actual.eigensystem.keys==reduced.keys
    @test actual.H.data ≈ Diagonal(collect(reduced.energies))
end

@testset "complete product partitions and synthesized empty factors" begin
    ax,bx,by=OperatorRef(:a,:x),OperatorRef(:b,:x),OperatorRef(:b,:y)
    A=QT.Qobj(ComplexF64[0 1 im;2 0 1;0 1 3])
    B=QT.sigmap();C=QT.sigmaz()
    repeated=((ax,false),(ax,true),(ax,false))
    rightx,righty=((bx,true),),((by,false),)
    left=CountedProducts(Dict{Tuple,QT.QuantumObject}(repeated=>A*A'*A),Dict{Tuple,Int}())
    right=CountedProducts(Dict{Tuple,QT.QuantumObject}(rightx=>B',righty=>C),Dict{Tuple,Int}())
    models=[(;owners=Set([:a]),operator_products=left),(;owners=Set([:b]),operator_products=right)]
    P=QuantumDevices.selection(3,[1,3])
    Q=QuantumDevices.projector_from_states([QT.basis(2,1)])
    Pdata,Qdata=copy(P.data),copy(Q.data)
    leftdata,rightdata=copy(left.data[repeated].data),copy(right.data[rightx].data)
    words=Set{Tuple}([(repeated...,rightx...),(repeated...,righty...),repeated,rightx,()])
    products=QuantumDevices.projected_products(models,[P,Q],words)
    @test left.reads==Dict(repeated=>1)
    @test right.reads==Dict(rightx=>1,righty=>1)
    @test !haskey(left.reads,()) && !haskey(right.reads,())
    @test products[()] == QT.tensor(QT.qeye(2),QT.qeye(1))
    @test products[repeated] ≈ QT.tensor(P'*(A*A'*A)*P,QT.qeye(1))
    @test products[(repeated...,rightx...)] ≈ QT.tensor(P'*(A*A'*A)*P,Q'*B'*Q)
    @test norm(products[repeated]-QT.tensor((P'*A*P)*(P'*A'*P)*(P'*A*P),QT.qeye(1)))>1e-5
    @test products[repeated].dims==([2,1],[2,1])
    @test P.data==Pdata && Q.data==Qdata
    @test left.data[repeated].data==leftdata
    @test right.data[rightx].data==rightdata
    # A purely empty request works even without any upstream product dictionary.
    empty_models=[(;owners=Set([:a]),operator_products=Dict())]
    @test QuantumDevices.projected_products(empty_models,[P],Set{Tuple}([()]))[()] == QT.qeye(2)
end
