using YAML

@testset "operator function IO and bounded legacy records" begin
    QD=QuantumDevices
    # Every builtin component declares portable identities for its own functions.
    for component in (make_qubit(:q), make_resonator(:r),
        make_transmon(:t; EC=.2, EJ=15.), make_flux_tunable_transmon(:f; EC=.2, EJmax=15.))
        for generator in values(component.operator_generators)
            record = QD._operator_function_record(generator)
            @test record !== nothing
            @test QD._restore_operator_function(record, "test") === generator
        end
    end
    K=Symbolics.variable(:portable_K)
    r=make_resonator(:r;dimension=4,frequency=5.)
    n=r.operators["n"]
    corrected=apply_correction(r,K/2*n*(n-1);parameters=Dict("K"=>(K,-.02)))
    device=QuantumDevice(:corrected);add_component!(device,corrected)
    add_recipe!(device,:leaf,Leaf(component=:r))
    folder=joinpath(mktempdir(),"corrected")
    save_device(folder,device)
    @test YAML.load_file(joinpath(folder,"device.yaml"))["schema"]==5
    restored=load_device(folder)
    @test restored.components[:r].dimension==4
    @test restored.components[:r].operator_generators["n"] === r.operator_generators["n"]
    @test isequal(restored.components[:r].hamiltonian,corrected.hamiltonian)
    @test build(restored,:leaf).spectrum ≈ build(device,:leaf).spectrum
    @test_throws ArgumentError QD._restore_operator_function(Dict("version"=>1,"id"=>"unknown"),"test")
    @test_throws ArgumentError QD._restore_operator_function(Dict("version"=>2,"id"=>"number"),"test")
    custom=make_resonator(:custom;dimension=2)
    custom.operator_generators["n"]=dimension->QT.num(dimension)
    custom_device=QuantumDevice(:custom);add_component!(custom_device,custom)
    @test_throws ArgumentError save_device(joinpath(mktempdir(),"custom"),custom_device)

    function legacy_context(d; recipes=true)
        local context_writer=QD.RecordWriter()
        local context_record=QD._encode_record(QD.capture_gate_context(d;include_recipes=recipes),context_writer,"context")
        local recipe_records=deepcopy(context_writer.recipes)
        legacy_generator_records!(context_record)
        foreach(legacy_generator_records!,values(recipe_records))
        return context_record,recipe_records
    end
    function decode_migrated(context,records,gates=Dict{String,Any}())
        QD._migrate_operator_records!(context,records,gates;schema=2)
        local migration_reader=QD.RecordReader(records,Dict(),Dict(),Set(),Dict(),Dict())
        return QD._decode_record(context,migration_reader,"context"),migration_reader
    end
    context,records=legacy_context(device)
    migrated,_=decode_migrated(context,records)
    @test migrated.components[:r].dimension==4
    @test migrated.recipes[:leaf].dimension==4
    @test migrated.components[:r].operator_generators["n"] === QD.number_operator
    @test build(QuantumDevice(migrated),:leaf).spectrum ≈ build(device,:leaf).spectrum

    # Legacy per-operator overrides could select a different implicit common default.
    override_device=QuantumDevice(:overrides)
    add_component!(override_device,make_resonator(:r;dimension=4))
    overrides=Dict(OperatorRef(:r,:a)=>QD.annihilation_operator,OperatorRef(:r,:n)=>QD.number_operator)
    add_recipe!(override_device,:leaf,Leaf(component=:r,generators=overrides))
    context,records=legacy_context(override_device)
    migrated,_=decode_migrated(context,records)
    @test migrated.components[:r].dimension==4
    @test migrated.recipes[:leaf].dimension==2

    # Partial overrides whose defaults differed were ambiguous, not component defaults.
    partial=QuantumDevice(:partial);add_component!(partial,make_resonator(:r;dimension=4))
    add_recipe!(partial,:leaf,Leaf(component=:r,generators=Dict(OperatorRef(:r,:n)=>QD.number_operator)))
    context,records=legacy_context(partial)
    @test_throws ArgumentError decode_migrated(context,records)
    add_recipe!(partial,:explicit,Leaf(component=:r,dimension=3,generators=Dict(OperatorRef(:r,:n)=>QD.number_operator)))
    delete!(partial.recipes,:leaf)
    context,records=legacy_context(partial)
    migrated,_=decode_migrated(context,records)
    @test migrated.recipes[:explicit].dimension==3

    # One interned recipe under two saved contexts must retain each old cutoff.
    small=QuantumDevice(:small);add_component!(small,make_resonator(:r;dimension=2))
    large=QuantumDevice(:large);add_component!(large,make_resonator(:r;dimension=5))
    writer=QD.RecordWriter()
    recipe=Leaf(component=:r)
    ref=QD._encode_record(recipe,writer,"recipe")
    contexts=[QD._encode_record(QD.capture_gate_context(d;include_recipes=false),writer,"context") for d in (small,large)]
    records=deepcopy(writer.recipes)
    foreach(legacy_generator_records!,contexts);foreach(legacy_generator_records!,values(records))
    # This named tuple has the same recipe/context structure as derived map provenance.
    derivations=Dict("small"=>Dict("tag"=>"named_tuple","keys"=>["recipe","context"],"values"=>Any[deepcopy(ref),contexts[1]]),
        "large"=>Dict("tag"=>"named_tuple","keys"=>["recipe","context"],"values"=>Any[deepcopy(ref),contexts[2]]))
    base,_=legacy_context(small;recipes=false)
    QD._migrate_operator_records!(base,records,derivations;schema=2)
    reader=QD.RecordReader(records,Dict(),Dict(),Set(),Dict(),Dict())
    a=QD._decode_record(derivations["small"],reader,"small")
    b=QD._decode_record(derivations["large"],reader,"large")
    @test a.recipe.dimension==2 && a.context.components[:r].dimension==2
    @test b.recipe.dimension==5 && b.context.components[:r].dimension==5
    @test derivations["small"]["values"][1]["id"] != derivations["large"]["values"][1]["id"]

    # A standalone historical artifact has no physical context of its own.
    # Its arrays do not justify borrowing the current device's implicit cutoff.
    context,records=legacy_context(small;recipes=false)
    artifact_writer=QD.RecordWriter()
    artifact=QD._encode_record(build(small,Leaf(component=:r)),artifact_writer,"artifact")
    foreach(legacy_generator_records!,values(artifact_writer.recipes))
    merge!(records,artifact_writer.recipes)
    @test_throws ArgumentError QD._migrate_operator_records!(context,records,Dict("model"=>artifact);schema=2)
    context,records=legacy_context(small;recipes=false)
    artifact_writer=QD.RecordWriter()
    artifact=QD._encode_record(build(small,Leaf(component=:r,dimension=2)),artifact_writer,"artifact")
    foreach(legacy_generator_records!,values(artifact_writer.recipes))
    merge!(records,artifact_writer.recipes)
    artifacts=Dict("model"=>artifact)
    QD._migrate_operator_records!(context,records,artifacts;schema=2)
    reader=QD.RecordReader(records,artifact_writer.arrays,Dict(),Set(),Dict(),Dict())
    @test QD._decode_record(artifacts["model"],reader,"model").recipe.definition.dimension==2

    # Known builtins ignored kwargs; custom configured functions cannot be guessed.
    context,records=legacy_context(device)
    leaf=only(values(records))
    leaf["fields"]["kwargs"]=Dict("tag"=>"named_tuple","keys"=>["ignored"],"values"=>[Dict("tag"=>"bool","value"=>true)])
    migrated,_=decode_migrated(context,records)
    @test migrated.recipes[:leaf].dimension==4
    context,records=legacy_context(override_device)
    leaf=only(values(records))
    leaf["fields"]["kwargs"]=Dict("tag"=>"named_tuple","keys"=>["scale"],"values"=>[QD._legacy_dimension_record(2)])
    first(leaf["fields"]["generators"]["entries"])[2]["fields"]["f"]=Dict("tag"=>"custom","id"=>"opaque")
    @test_throws ArgumentError decode_migrated(context,records)
    leaf=Dict("kind"=>"Leaf","fields"=>Dict("component"=>Dict("tag"=>"symbol","value"=>"r"),
        "dimension"=>Dict("tag"=>"nothing"),"generators"=>Dict("tag"=>"dict","entries"=>Any[])))
    @test_throws ArgumentError QD._migrate_leaf_generators!(leaf,nothing,"missing")
end
