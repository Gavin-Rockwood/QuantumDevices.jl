@testset "registered recipe snapshots and editable drafts" begin
    device = QuantumDevice(:snapshots)
    q = add_component!(device, make_qubit(:q; frequency=5.0))
    child = Recipe(device, [q]; states=Dict(q=>[0,1]), parameters=Dict(q.variables["ν"]=>5.2))
    coupling = Symbolics.variable(:snapshot_shift)
    original_parameters=Dict(coupling=>.1)
    interaction = Interaction(:shift, (coupling)*child.inputs.q.operators.z)
    merge!(child.overrides, original_parameters)
    add_interaction!(child, interaction)
    saved_child = add_recipe!(device, :child, child)
    target = Recipe(device, [q]; inputs=(child=:child,), states=Dict(:child=>[1,0]))
    saved_target = add_recipe!(device, :target, target)
    reference = build(device, :target)
    child.overrides[q.variables["ν"]] = 9.0
    child.states[:q][1] = 1
    empty!(child.interactions)
    empty!(original_parameters)
    child.overrides[coupling] = .9
    @test_throws Exception (child.dressing_options.nsteps = 9)
    @test saved_child.overrides[q.variables["ν"]] == 5.2
    @test collect(saved_child.states[:q]) == [0,1]
    @test length(saved_child.interactions) == 1
    @test saved_child.overrides[coupling] == .1
    @test saved_child.dressing_options.nsteps == 20
    @test build(device,:target).spectrum ≈ reference.spectrum
    @test saved_target !== target
    @test saved_target.inputs.child.recipe === saved_child ||
        get_hamiltonian(device,saved_target.inputs.child.recipe) == get_hamiltonian(device,saved_child)
    @test_throws Exception (saved_child.overrides[q.variables["ν"]] = 8.0)
    @test_throws Exception (saved_child.states[:q][1] = 1)
    @test_throws Exception push!(saved_child.interactions,interaction)
    @test_throws Exception (saved_child.overrides[coupling] = .2)
    @test_throws Exception (saved_child.dressing_options.nsteps = 8)
    trial = draft(saved_target)
    trial.overrides[q.variables["ν"]] = 6.0
    @test trial.overrides isa Dict
    @test trial.inputs.child.recipe.overrides isa Dict
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
    @test reference.recipe.definition.inputs.child.recipe.overrides[q.variables["ν"]] == 5.2
    model = build(device,trial)
    old_spectrum = copy(model.spectrum)
    trial.overrides[q.variables["ν"]] = 8.0
    trial.states[:child][1] = 0
    @test model.resolved_parameters[q.variables["ν"]]==6.
    @test model.recipe.definition.overrides[q.variables["ν"]]==6.
    @test model.spectrum==old_spectrum
    @test_throws Exception (model.resolved_parameters[q.variables["ν"]]=99.)
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
    generators=Dict(OperatorRef(:q,:x)=>q.operator_generators["x"])
    leaf=Leaf(component=:q,generators=generators)
    saved=add_recipe!(d,:leaf,leaf)
    empty!(generators)
    @test length(saved.generators)==1
    @test_throws Exception empty!(saved.generators)
    trial=draft(saved)
    empty!(trial.generators)
    @test length(saved.generators)==1
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

@testset "generator configuration belongs to ordinary closures" begin
    d = QuantumDevice(:closure_configuration)
    q = add_component!(d, make_qubit(:q))
    payload = [1.0, 2.0]
    generator = dimension -> payload[1] * QT.sigmax()
    recipe = Leaf(component=:q, generators=Dict(OperatorRef(:q,:x)=>generator))
    saved = add_recipe!(d,:saved,recipe)
    @test saved.generators[OperatorRef(:q,:x)] === generator
    @test build(d,:saved).operators[OperatorRef(:q,:x)] ≈ QT.sigmax()
    # Snapshot ownership does not pretend to freeze externally captured state.
    payload[1] = 2.0
    @test build(d,:saved).operators[OperatorRef(:q,:x)] ≈ 2 * QT.sigmax()
end
