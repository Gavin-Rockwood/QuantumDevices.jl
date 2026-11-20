@testset "dressed data energy order, tensor keys and explicit state ownership" begin
    d=QuantumDevice(:tracked_models)
    q1=add_component!(d,make_resonator(:q1;dimension=3,frequency=2.1))
    q2=add_component!(d,make_resonator(:q2;dimension=3,frequency=3.7))
    q3=add_component!(d,make_resonator(:q3;dimension=2,frequency=5.3))
    leaf=build(d,Leaf(component=:q1))
    @test leaf.eigensystem.keys==((0,),(1,),(2,))
    @test leaf.eigensystem.components==(:q1,)
    @test leaf.energies(1)≈2.1
    @test (@inferred leaf.energies(1))≈2.1
    @test (@inferred leaf.states(1))≈input_basis_map(leaf)*QT.basis(3,1)
    recipe=Recipe(d,[q1,q2,q3];states=Dict(q1=>[2,0],q2=>[1,0],q3=>0:1))
    model=build(d,recipe)
    @test model.eigensystem.components==(:q1,:q2,:q3)
    @test length(model.eigensystem.keys)==8
    @test model.states.tracking===model.energies.tracking===model.eigensystem.tracking
    @test getfield(model.states, :_storage) isa Vector{<:QT.QuantumObject}
    @test !hasproperty(model, :rotation)
    @test :rotation ∉ fieldnames(typeof(model))
    @test !isdefined(QuantumDevices, :RotationStates)
    @test model.eigensystem.tracking.saved_steps==(2,)
    @test model.eigensystem.tracking.total_steps==2 # Exactly zero V skips repeated ramp work.
    @test model.eigensystem.tracking.indices==collect(0:7)
    @test :states in propertynames(model) && :energies in propertynames(model)
    @test model.spectrum===model.energies || model.spectrum==model.energies
    @test model.labeling===model.eigensystem.tracking
    for i in 0:1,j in 0:1,k in 0:1
        key=([2,0][i+1],[1,0][j+1],k)
        expected=QT.tensor(QT.basis(2,i),QT.basis(2,j),QT.basis(2,k))
        @test abs2(dot(expected,model.states(key)))≈1
        @test model.states(key) === state(model,key)
        @test stateenergy(model,key)≈2.1*key[1]+3.7*key[2]+5.3*k
        @test stateindex(model,key)==findfirst(==(key),model.eigensystem.keys)-1
        @test model.energies(key)==model.energies(key...)
    end
    @test state(model,2,1,0)≈state(model,(2,1,0))
    @test_throws KeyError model.states(2,1)
    @test_throws KeyError model.states(1,0,0)
    @test states(model)[1] === model.states[1]
    pair=Recipe(d,[q2,q1];states=Dict(q1=>[2,0],q2=>[1,0]))
    pairmodel=build(d,pair)
    nested=build(d,Recipe(d,[q3,q1,q2];inputs=(group=pair,),states=Dict(:group=>[2,0],q3=>[1])))
    @test nested.eigensystem.components==(:q2,:q1,:q3)
    @test Set(nested.eigensystem.keys)==Set((key...,1) for key in pairmodel.eigensystem.keys[[3,1]])
    for key in nested.eigensystem.keys
        @test nested.energies(key)≈3.7*key[1]+2.1*key[2]+5.3*key[3]
    end
    off=Recipe(d,[q1,q2];dressing=false)
    unlabeled=build(d,off)
    @test isempty(unlabeled.eigensystem.keys) && isnothing(unlabeled.labeling)
    @test isnothing(unlabeled.eigensystem.tracking.accuracy)
    @test isnothing(unlabeled.eigensystem.tracking.diagnostics)
    @test unlabeled.energies[1]≈0
    @test unlabeled.states[1]≈input_basis_map(unlabeled)*QT.basis(9,0)
    @test_throws ArgumentError unlabeled.energies(0,0)
    @test_throws ArgumentError state(unlabeled,0,0)
    @test_throws ArgumentError build(d,Recipe(d,[q1,q2];inputs=(group=off,)))
    @test isnothing(build(d,Recipe(d,[q1,q2];inputs=(group=off,),dressing=false)).labeling)
    @test_throws ArgumentError DressedStates(leaf.eigensystem.tracking,(:q1,:q2))
    @test_throws ArgumentError DressedStates(leaf.eigensystem.tracking,(:q1,:q1))
    badkeys=TrackingResult(((0,),(true,)),(1,),1,[QT.basis(2,0),QT.basis(2,1)],(energies=[0.,1.],),[0,1],
        (overlaps=ones(2),),[LabelDiagnostic(1.,0.,false,false),LabelDiagnostic(1.,0.,false,false)])
    @test_throws ArgumentError DressedStates(badkeys,(:q1,))
    @test occursin("8 states",repr(model.eigensystem))
    @test occursin("physical component order",repr(MIME"text/plain"(),model.eigensystem))
    @test !occursin("Dict",repr(MIME"text/plain"(),model.eigensystem))
    legacy=NumericalModel(model.recipe,model.H,QuantumDevices.untracked_eigensystem(collect(model.spectrum),input_basis_map(model),model.eigensystem.components),model.operators,model.hamiltonian,model.owners,model.operator_products)
    @test legacy.spectrum==model.spectrum && isnothing(legacy.labeling)
    @test !hasproperty(legacy, :rotation)
    @test all(isapprox.(legacy.states, model.states))
end

@testset "ramp continuation, inherited accuracy and phase-consistent eigenvectors" begin
    d=QuantumDevice(:avoided_crossing)
    q=add_component!(d,make_qubit(:q;frequency=2.))
    function crossing_recipe(options=(;); dressing=true)
        fixture_recipe=Recipe(d,[q];dressing, dressing_options=options)
        add_interaction!(fixture_recipe,Interaction(:ramp,-2*fixture_recipe.inputs.q.operators.z+.1*fixture_recipe.inputs.q.operators.x))
        fixture_recipe
    end
    r=crossing_recipe((;nsteps=101,trajectory=identity))
    model=build(d,r)
    finer=build(d,crossing_recipe((;nsteps=201,trajectory=identity)))
    endpoint=build(d,crossing_recipe(;dressing=false))
    @test stateindex(model,0)==0
    @test model.eigensystem.keys==finer.eigensystem.keys
    @test abs2(dot(QT.basis(2,0),model.states[2]))>.9
    @test model.spectrum≈endpoint.spectrum
    @test make_hamiltonian(model) ≈ make_hamiltonian(endpoint)
    @test input_basis_map(model)'*input_basis_map(model)≈QT.qeye(2)
    hfinal=-QuantumDevices.spectrum_operator([-1.,1.])+.1*QT.sigmax()
    @test hfinal*input_basis_map(model)≈input_basis_map(model)*QuantumDevices.spectrum_operator(model.spectrum)
    for i in 1:2
        @test model.states[i]≈input_basis_map(model)*QT.basis(2,i-1)
        @test hfinal*model.states[i]≈model.energies[i]*model.states[i]
    end
    @test model.eigensystem.tracking.total_steps==101
    @test model.eigensystem.tracking.saved_steps==(101,)
    combined=build(d,Recipe(d,[q];inputs=(child=r,)))
    for key in model.eigensystem.keys
        before=model.eigensystem.tracking.diagnostics(key)
        after=combined.eigensystem.tracking.diagnostics(key)
        @test after.minimum_overlap<=before.minimum_overlap+1e-12
        @test after.cumulative_log_overlap≈before.cumulative_log_overlap
        @test combined.eigensystem.tracking.accuracy.minimum_overlaps(key)==after.minimum_overlap
        @test combined.eigensystem.tracking.accuracy.cumulative_log_overlaps(key)==after.cumulative_log_overlap
    end
    degenerateq=add_component!(d,make_qubit(:degenerate;frequency=0.))
    degenerate=build(d,Recipe(d,[degenerateq]))
    @test all(x->x.degenerate,degenerate.eigensystem.tracking.diagnostics)
end

