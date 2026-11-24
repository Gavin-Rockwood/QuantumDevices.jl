using YAML, SHA
@testset "Complete device folders" begin
    QD=QuantumDevices
    device=QuantumDevice(:portable)
    legacy=QuantumDevice(device.name,device.components,device.interactions,device.recipes,device.metadata)
    @test isempty(legacy.gates)
    q=make_resonator(:q;dimension=2,frequency=1.)
    add_component!(device,q)
    recipe=Recipe(device,[q];initial_dims=Dict(q=>2),states=Dict(q=>0:1))
    add_recipe!(device,:model,recipe)
    gate=QuantumGate(:X,GateRecipe(recipe;controls=[Control(:drive,q.variables["ν"],Pulse(RampedFlattop(.01;duration=.1);duration=.1,amplitude=.1,offset=1.))]))
    @test add_gate!(device,gate)===gate
    @test_throws ArgumentError add_gate!(device,gate)
    @test move_gate!(device,:X,(:X,:small))===gate
    @test get_gate(device,(:X,:small))===gate
    @test move_gate!(device,(:X,:small),(:X,:small))===gate
    @test_throws ArgumentError move_gate!(device,(:X,:small),:X)
    @test get_gate(device,(:X,:small))===gate
    @test_throws ArgumentError get_gate(device,:X)
    @test_throws ArgumentError replace_gate!(device,:absent,gate)
    replacement=QuantumGate(:DifferentName,gate.recipe)
    @test replace_gate!(device,(:X,:small),replacement)===replacement
    @test get_gate(device,(:X,:small)).name==:DifferentName
    add_gate!(device,:alias,replacement)
    @test_throws ArgumentError attach_artifact!(replacement,:bad,()->0)
    root=mktempdir();folder=joinpath(root,"device")
    @test save_device(folder,device)==folder
    @test isfile(joinpath(folder,"device.yaml"))
    loaded=load_device(folder)
    @test loaded.name==device.name
    @test get_gate(loaded,:alias)===get_gate(loaded,(:X,:small))
    @test isequal(loaded.parameters,device.parameters)
    @test_throws ArgumentError save_device(folder,device)
    before=read(joinpath(folder,"device.yaml"))
    @test save_device(folder,device;overwrite=true)==folder
    @test_throws ArgumentError QD._bundle_path(folder,"../outside")
    @test_throws ArgumentError QD._bundle_path(folder,"/outside")
    @test_throws ArgumentError QD._bundle_path(folder,"a/./b")
    # Equivalent captured/draft recipes collapse, including refreshed scope identities.
    w=QD.RecordWriter()
    first=QD.intern_recipe!(w,recipe,"first")
    second=QD.intern_recipe!(w,QD.snapshot_recipe(device,recipe),"second")
    @test first["id"]==second["id"]
    changed=setpath(recipe,"max_dimension",999)
    @test QD.intern_recipe!(w,changed,"changed")["id"]!=first["id"]
    # Canonical maps ignore insertion order, but exact values/configuration remain semantic.
    left=Leaf(component=:q,dimension=2,parameters=Dict(q.variables["ν"]=>1.))
    right=Leaf(component=:q,dimension=2,parameters=Dict(q.variables["ν"]=>1.))
    @test QD.intern_recipe!(w,left,"left")["id"]==QD.intern_recipe!(w,right,"right")["id"]
    typed=Leaf(component=:q,dimension=2,parameters=Dict(q.variables["ν"]=>Int32(1)))
    @test QD.intern_recipe!(w,typed,"typed")["id"]!=QD.intern_recipe!(w,left,"left")["id"]
    defaults1=Leaf(component=:q,dimension=2,parameters=Dict(q.variables["ν"]=>1.))
    defaults2=Leaf(component=:q,dimension=2,parameters=Dict(q.variables["ν"]=>2.))
    @test QD.intern_recipe!(w,defaults1,"defaults1")["id"]!=QD.intern_recipe!(w,defaults2,"defaults2")["id"]
    collision=QD.RecordWriter();id=QD.intern_recipe!(collision,recipe,"original")["id"]
    collision.recipes[id]=Dict("different"=>"canonical bytes")
    @test_throws ArgumentError QD.intern_recipe!(collision,recipe,"collision")
    # Physical defaults and scalar ASTs survive without generating operators on load.
    physical=QuantumDevice(:physical)
    transmon=make_flux_tunable_transmon(:t,.2,12.;dimension=7,flux=.1)
    add_component!(physical,transmon)
    save_device(joinpath(root,"physical"),physical)
    restored=load_device(joinpath(root,"physical"))
    @test isequal(restored.components[:t].hamiltonian,transmon.hamiltonian)
    @test isequal(restored.parameters,physical.parameters)
    # Unassigned registered interaction variables survive portable declarations.
    unresolved = Symbolics.variable(:unresolved)
    add_interaction!(physical,Interaction(:unset,unresolved*transmon.operators["charge"]))
    save_device(joinpath(root,"unresolved"),physical)
    restored_unresolved = load_device(joinpath(root,"unresolved"))
    @test isequal(restored_unresolved.variables,physical.variables)
    @test isequal(restored_unresolved.parameters,physical.parameters)
    @test isempty(restored_unresolved.parameters["unresolved"])
    # A tiny calibration owns complete reference arrays and can replay/continue portably.
    prepared=build(device,replacement)
    A=Symbolics.variable(:portable_A)
    problem=CalibrationProblem(prepared,Matrix{ComplexF64}(I,2,2);variables=[A=>.1],links=Dict("drive/amplitude"=>A),
        states=[(0,),(1,)],lb=[0.],ub=[.3])
    native=(;u=[.1],retcode=QD.SciMLBase.ReturnCode.MaxIters,objective=.2)
    result=CalibrationResult(problem,native)
    calibrated=QuantumGate(:historical,result;allow_unsuccessful=true)
    add_gate!(device,:historical,calibrated)
    evolution=evolve(result)
    attach_artifact!(calibrated,:evolution,evolution)
    attach_artifact!(calibrated,:model,prepared.model)
    @test_throws ArgumentError attach_artifact!(calibrated,:model,prepared.model)
    native.u[1]=.2 # Summary u and authoritative fitted_values deliberately differ.
    save_device(folder,device;overwrite=true,native_checkpoint=true)
    loaded=load_device(folder);saved=get_gate(loaded,:historical)
    @test saved.calibration.fitted_values==(.1,)
    @test saved.calibration.solution.u==[.2]
    @test saved.calibration.solution.retcode==QD.SciMLBase.ReturnCode.MaxIters
    @test saved.calibration.provenance.prepared.model.H.data==prepared.model.H.data
    @test saved.calibration.provenance.prepared.model.operators==prepared.model.operators
    @test saved.artifacts[:model].H.dimensions==prepared.model.H.dimensions
    @test gate_matrix(evolve(saved.calibration)).data≈gate_matrix(evolution).data
    @test saved.artifacts[:evolution].propagators[end].data≈evolution.propagators[end].data
    continuation=CalibrationProblem(saved.calibration)
    @test continuation.initial_values==(.1,)
    @test_throws ArgumentError load_device(folder;load_native=true)
    native_loaded=load_device(folder;load_native=true,trusted_native=true)
    @test get_gate(native_loaded,:historical).calibration.solution isa typeof(native)
    @test get_gate(native_loaded,:historical).calibration.fitted_values==(.1,)
    actual=solve(problem;maxiters=2)
    add_gate!(device,:native,QuantumGate(:native,actual;allow_unsuccessful=true))
    save_device(folder,device;overwrite=true,native_checkpoint=true)
    @test get_gate(load_device(folder;load_native=true,trusted_native=true),:native).calibration.solution isa typeof(actual.solution)
    @test get_gate(load_device(folder),:native).calibration.solution isa NamedTuple
    checkpoint=joinpath(folder,"checkpoints/native.jld2")
    checkpoint_bytes=read(checkpoint);rm(checkpoint)
    @test get_gate(load_device(folder),:native).calibration.fitted_values==actual.fitted_values
    @test_throws ArgumentError load_device(folder;load_native=true,trusted_native=true)
    write(checkpoint,checkpoint_bytes)
    # Independent working and historical owners despite shared definition IDs.
    saved.recipe=setpath(saved.recipe,"drive/amplitude",.25)
    @test getpath(saved.calibration.recipe,"drive/amplitude")==.1
    @test getpath(get_gate(loaded,:alias).recipe,"drive/amplitude")==.1
    # Portable numerical records do not invoke arbitrary Julia reconstruction.
    bad=QuantumDevice(:bad;metadata=(;callback=()->0))
    @test_throws ArgumentError save_device(folder,bad;overwrite=true)
    @test isfile(joinpath(folder,"device.yaml"))
    # Fail before returning a partial device when a declared file is corrupted.
    manifest=YAML.load_file(joinpath(folder,"device.yaml"))
    payload=joinpath(folder,manifest["payload"])
    open(payload,"a") do io;write(io,UInt8(1));end
    @test_throws ArgumentError load_device(folder)
end
