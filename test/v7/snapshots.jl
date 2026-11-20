@testset "registered recipe snapshots and editable drafts" begin
    device = QuantumDevice(:snapshots)
    q = add_component!(device, make_qubit(:q; frequency=5.0))
    child = Recipe(device, [q]; states=Dict(q=>[0,1]), parameters=Dict(q.variables["ν"]=>5.2))
    coupling = Symbolics.variable(:snapshot_shift)
    original_parameters=Dict(coupling=>.1)
    interaction = Interaction(:shift, (coupling)*child.inputs.q.operators.z; parameters=original_parameters)
    add_interaction!(child, interaction)
    saved_child = add_recipe!(device, :child, child)
    target = Recipe(device, [q]; inputs=(child=:child,), states=Dict(:child=>[1,0]))
    saved_target = add_recipe!(device, :target, target)
    reference = build(device, :target)
    child.parameters[q.variables["ν"]] = 9.0
    child.states[:q][1] = 1
    empty!(child.interactions)
    empty!(interaction.parameters)
    @test_throws Exception (child.dressing_options.nsteps = 9)
    @test saved_child.parameters[q.variables["ν"]] == 5.2
    @test collect(saved_child.states[:q]) == [0,1]
    @test length(saved_child.interactions) == 1
    @test length(saved_child.interactions[1].parameters) == 1
    @test saved_child.dressing_options.nsteps == 20
    @test build(device,:target).spectrum ≈ reference.spectrum
    @test saved_target !== target
    @test saved_target.inputs.child.recipe === saved_child ||
        get_hamiltonian(device,saved_target.inputs.child.recipe) == get_hamiltonian(device,saved_child)
    @test_throws Exception (saved_child.parameters[q.variables["ν"]] = 8.0)
    @test_throws Exception (saved_child.states[:q][1] = 1)
    @test_throws Exception push!(saved_child.interactions,interaction)
    @test_throws Exception empty!(saved_child.interactions[1].parameters)
    @test_throws Exception (saved_child.dressing_options.nsteps = 8)
    trial = draft(saved_target)
    trial.parameters[q.variables["ν"]] = 6.0
    @test trial.parameters isa Dict
    @test trial.inputs.child.recipe.parameters isa Dict
    unchanged = draft(saved_target)
    @test build(device,unchanged).spectrum ≈ reference.spectrum
    @test repr(MIME"text/plain"(),graph(device,unchanged)) == repr(MIME"text/plain"(),graph(device,saved_target))
    local_trial = draft(saved_child)
    @test local_trial.scope != saved_child.scope
    @test only(QuantumDevices.refs(local_trial.interactions[1].expression,:op)).scope == local_trial.scope
    add_interaction!(local_trial,Interaction(:new_shift,.2*local_trial.inputs.q.operators.z))
    @test length(local_trial.interactions)==2 && length(saved_child.interactions)==1
    replacement = replace_recipe!(device,:child,Recipe(device,[q];parameters=Dict(q.variables["ν"]=>7.)))
    @test replacement === device.recipes[:child]
    @test build(device,:target).spectrum ≈ reference.spectrum
    @test reference.recipe.definition.inputs.child.recipe.parameters[q.variables["ν"]] == 5.2
    model = build(device,trial)
    old_spectrum = copy(model.spectrum)
    trial.parameters[q.variables["ν"]] = 8.0
    trial.states[:child][1] = 0
    @test model.parameters[q.variables["ν"]]==6.
    @test model.recipe.definition.parameters[q.variables["ν"]]==6.
    @test model.spectrum==old_spectrum
    @test_throws Exception (model.parameters[q.variables["ν"]]=99.)
    @test_throws Exception push!(model.recipe.inputs,:oops=>model.recipe)
    @test_throws ArgumentError add_recipe!(device,:target,Recipe(device,[q]))
    bad = Combine(inputs=[Input(:missing,:missing,[0])])
    @test_throws ArgumentError add_recipe!(device,:failed,bad)
    @test !haskey(device.recipes,:failed)
    @test_throws ArgumentError replace_recipe!(device,:target,bad)
    @test device.recipes[:target] === saved_target
    @test_throws ArgumentError replace_recipe!(device,:absent,unchanged)
end

@testset "leaf options and legacy interaction capture" begin
    d=QuantumDevice(:legacy_snapshots)
    q=add_component!(d,make_qubit(:q))
    matrix = [1. 2.; 3. 4.]
    kwargs=(nested=Dict(:levels=>[0,1]), matrix=matrix,)
    generators=Dict(OperatorRef(:q,:x)=>q.operator_generators["x"])
    leaf=Leaf(component=:q,kwargs=kwargs,generators=generators)
    saved=add_recipe!(d,:leaf,leaf)
    kwargs.nested[:levels][1]=9
    matrix[1,1]=99.
    empty!(generators)
    @test saved.kwargs.nested[:levels][1]==0
    @test saved.kwargs.matrix[1,1]==1.
    @test length(saved.generators)==1
    @test_throws Exception (saved.kwargs.nested[:levels][1]=3)
    trial=draft(saved)
    @test trial.kwargs.matrix isa Matrix{Float64}
    trial.kwargs.nested[:levels][1]=2
    @test saved.kwargs.nested[:levels][1]==0
    add_interaction!(d,Interaction(:physical,.1*op(:q,:z)))
    combine=Combine(inputs=[Input(:q,:leaf,[0,1])],interactions=[:physical])
    frozen=add_recipe!(d,:combine,combine)
    before=build(d,:combine).spectrum
    d.interactions[:physical]=Interaction(:physical,.9*op(:q,:z))
    @test build(d,:combine).spectrum≈before
    @test build(d,draft(frozen)).spectrum≈before
    combine.inputs[1].states[1]=1
    @test collect(frozen.inputs[1].states)==[0,1]
end

@testset "generator runtime keyword copies" begin
    d = QuantumDevice(:keyword_copy)
    q = add_component!(d, make_qubit(:q))
    values = [1.0, 2.0]
    generator = OperatorGenerator((dimension; payload::Vector{Float64}) -> begin
        payload[1] = 99.0
        QT.sigmax()
    end, 2)
    recipe = Leaf(component=:q, generators=Dict(OperatorRef(:q,:x)=>generator), kwargs=(payload=values,))
    saved = add_recipe!(d,:saved,recipe)
    @test build(d,:saved).operators[OperatorRef(:q,:x)] isa QT.QuantumObject
    @test values==[1.,2.]
    @test collect(saved.kwargs.payload)==[1.,2.]
    @test build(d,:saved).recipe.definition.kwargs.payload[1]==1.
end
