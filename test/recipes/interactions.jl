import JLD2

@testset "physical interactions keyword" begin
    d = QuantumDevice(:physical_interactions)
    r = add_component!(d,make_resonator(:r;dimension=4,frequency=3.))
    q = add_component!(d,make_qubit(:q;frequency=5.))
    g = Symbolics.variable(:physical_g)
    a = r.operators["a"]
    source = Interaction(:product,g*(a*a')*q.operators["z"];
        metadata=(notes=["original"],))
    automatic = Recipe(d,[r,:q];interactions=[source],parameters=Dict(g=>.07),dressing=false)
    manual = Recipe(d,[r,q];parameters=Dict(g=>.07),dressing=false)
    ma = manual.inputs.r.operators.a
    add_interaction!(manual,Interaction(:product,g*(ma*ma')*manual.inputs.q.operators.z))
    @test build(d,automatic).spectrum ≈ build(d,manual).spectrum
    @test all(ref -> ref isa OperatorRef,QuantumDevices.refs(source.expression,:op))
    @test only(automatic.interactions) !== source
    @test only(automatic.interactions).variables !== source.variables
    @test !hasproperty(only(automatic.interactions), :parameters)
    push!(only(automatic.interactions).metadata.notes,"captured")
    @test source.metadata.notes == ["original"]
    # Registered device interactions are selected explicitly and counted once.
    add_interaction!(d,source)
    @test build(d,Recipe(d,[r,q];interactions=[source],parameters=Dict(g=>.07),dressing=false)).spectrum ≈ build(d,manual).spectrum
    @test isempty(Recipe(d,[r,q]).interactions)
    group = Recipe(d,[r,q];states=Dict(r=>[0,1],q=>[0,1]),dressing=false)
    grouped = Recipe(d,[r,q];inputs=(pair=group,),interactions=[source],parameters=Dict(g=>.07),dressing=false)
    grouped_manual = Recipe(d,[r,q];inputs=(pair=group,),parameters=Dict(g=>.07),dressing=false)
    ga = grouped_manual.inputs.pair.operators.r.a
    add_interaction!(grouped_manual,Interaction(:product,g*(ga*ga')*grouped_manual.inputs.pair.operators.q.z))
    @test build(d,grouped).spectrum ≈ build(d,grouped_manual).spectrum
    @test all(ref -> ref.input == :pair,QuantumDevices.refs(only(grouped.interactions).expression,:op))
    @test_throws ArgumentError Recipe(d,[r];interactions=[source])
    @test_throws ArgumentError Recipe(d,[r,q];interactions=[Interaction(:bad,op(:r,:missing))])
    @test_throws ArgumentError Recipe(d,[r,q];interactions=[source,source])
    @test_throws ArgumentError Recipe(d,[r,q];inputs=(first=group,second=group),interactions=[source])
    @test_throws ArgumentError Recipe(d,[r,q];interactions=[:product])
    @test_throws ArgumentError Recipe(d,[r,q];interactions=[only(manual.interactions)])
    missing = Symbolics.variable(:physical_missing)
    term = Interaction(:missing,missing*q.operators["z"])
    unresolved = Recipe(d,[q];interactions=[term],dressing=false)
    @test_throws ArgumentError build(d,unresolved)
    @test build(d,unresolved;parameters=Dict(missing=>.2)).resolved_parameters[missing] == .2
    @test build(d,Recipe(d,[q];interactions=[term],parameters=Dict(missing=>.3))).resolved_parameters[missing] == .3
    assigned = Recipe(d,[q];interactions=[Interaction(:assigned,g*q.operators["z"])],parameters=Dict(g=>9.))
    @test build(d,assigned).resolved_parameters[g] == 9.
    frequency = q.variables["ν"]
    component_override = Recipe(d,[q];interactions=[Interaction(:component_shift,
        frequency*q.operators["z"])],parameters=Dict(frequency=>9.))
    @test build(d,component_override).resolved_parameters[frequency] == 9.
    # A parent assignment overrides the stored upstream operating point.
    upstream = Recipe(d,[q];interactions=[Interaction(:upstream,missing*q.operators["z"])],parameters=Dict(missing=>.1))
    overridden = Recipe(d,[q];inputs=(child=upstream,),
        interactions=[Interaction(:parent,missing*q.operators["z"])],parameters=Dict(missing=>9.))
    @test build(d,overridden).resolved_parameters[missing] == 9.
    add_recipe!(d,:saved_physical,grouped)
    saved = d.recipes[:saved_physical]
    editable = draft(saved)
    @test build(d,saved).spectrum ≈ build(d,grouped).spectrum
    @test build(d,editable).spectrum ≈ build(d,grouped).spectrum
    @test editable.scope != grouped.scope
    @test all(ref -> ref.scope == editable.scope,
        QuantumDevices.refs(only(editable.interactions).expression,:op))
    mktempdir() do directory
        path = joinpath(directory,"recipe.jld2")
        JLD2.jldsave(path;saved)
        restored = JLD2.load(path,"saved")
        @test build(d,restored).spectrum ≈ build(d,grouped).spectrum
        @test build(d,draft(restored)).spectrum ≈ build(d,grouped).spectrum
    end
end

@testset "automatic registered interactions" begin
    d = QuantumDevice(:auto_interactions)
    q1, coupler, q2, outside = [add_component!(d, make_qubit(name; frequency=f))
        for (name,f) in zip((:q1,:coupler,:q2,:outside),(4.,6.,5.,7.))]
    x1, xc, x2, xo = [c.operators["x"] for c in (q1,coupler,q2,outside)]
    # Intentionally register in a different order from the captured order.
    b = add_interaction!(d, Interaction(:b, .02*xc*x2))
    a = add_interaction!(d, Interaction(:a, .03*x1*xc))
    add_interaction!(d, Interaction(:outside_coupling, .01*x1*xo))
    hyper = add_interaction!(d, Interaction(:hyper, .04*x1*xc*x2))
    add_interaction!(d, Interaction(:single, .01*x1*x1))
    add_interaction!(d, Interaction(:scalar, QuantumDevices.val(0.1)))
    pair = Recipe(d, [q1,:coupler]; interactions=:auto, dressing=false)
    @test getfield.(pair.interactions,:name) == [:a]
    allthree = Recipe(d, [q1,coupler,q2]; interactions=:auto, dressing=false)
    @test getfield.(allthree.interactions,:name) == [:a,:b,:hyper]
    @test isempty(Recipe(d, [q1,q2]; interactions=:auto).interactions)
    @test isempty(Recipe(d, [q1,coupler,q2]).interactions)
    @test all(ref -> ref.scope == allthree.scope,
        Iterators.flatten(QuantumDevices.refs(i.expression,:op) for i in allthree.interactions))
    manual = Recipe(d, [q1,coupler,q2]; interactions=[a,b,hyper], dressing=false)
    @test build(d,allthree).spectrum ≈ build(d,manual).spectrum
    grouped = Recipe(d, [q1,coupler,q2]; inputs=(pair=pair,), interactions=:auto, dressing=false)
    @test getfield.(grouped.interactions,:name) == [:b,:hyper]
    grouped_manual = Recipe(d, [q1,coupler,q2]; inputs=(pair=pair,), interactions=[b,hyper], dressing=false)
    @test build(d,grouped).spectrum ≈ build(d,grouped_manual).spectrum
    @test build(d,grouped).spectrum ≈ build(d,manual).spectrum
    @test Set(ref.input for ref in QuantumDevices.refs(last(grouped.interactions).expression,:op)) == Set([:pair,:q2])
    @test isempty(Recipe(d, [q1,coupler]; inputs=(pair=pair,), interactions=:auto).interactions)
    # Even an uncoupled grouped input owns its internal terms.
    bare_pair = Recipe(d, [q1,coupler]; dressing=false)
    @test isempty(Recipe(d, [q1,coupler]; inputs=(pair=bare_pair,), interactions=:auto).interactions)
    add_interaction!(d, Interaction(:later, .01*x1*x2))
    @test getfield.(allthree.interactions,:name) == [:a,:b,:hyper]
    @test getfield.(Recipe(d,[q1,coupler,q2];interactions=:auto).interactions,:name) == [:a,:b,:hyper,:later]
    saved = add_recipe!(d,:auto_saved,grouped)
    @test build(d,saved).spectrum ≈ build(d,grouped).spectrum
    editable = draft(saved)
    @test getfield.(editable.interactions,:name) == [:b,:hyper]
    @test build(d,editable).spectrum ≈ build(d,grouped).spectrum
    @test all(ref -> ref.scope == editable.scope,
        Iterators.flatten(QuantumDevices.refs(i.expression,:op) for i in editable.interactions))
    mktempdir() do directory
        path = joinpath(directory,"auto_recipe.jld2")
        JLD2.jldsave(path; saved)
        restored = JLD2.load(path,"saved")
        @test getfield.(restored.interactions,:name) == [:b,:hyper]
        @test build(d,draft(restored)).spectrum ≈ build(d,grouped).spectrum
    end
    @test_throws ArgumentError Recipe(d,[q1,coupler];interactions=:unknown)
end
