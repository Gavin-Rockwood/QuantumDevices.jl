@testset "leaf input capture persistence and contextual legacy migration" begin
    QD=QuantumDevices
    d=QuantumDevice(:input_capture_io)
    q=add_component!(d,make_qubit(:q;frequency=5.))
    ν=q.variables["ν"]
    add_recipe!(d,:base,Recipe(d,[q]))
    q.parameters["ν"]=9.
    folder=joinpath(mktempdir(),"captured")
    save_device(folder,d)
    restored=load_device(folder)
    @test restored.components[:q].parameters["ν"]==9.
    @test haskey(restored.recipes[:base].parameters[string(ν)],5.)
    @test build(restored,:base).resolved_parameters[ν]==5.

    # Verified old records lack the leaf field. Reused recipe IDs are bound
    # independently under the current and historical component contexts.
    writer=QD.RecordWriter()
    definition=QD._encode_record(QD.capture_gate_context(d),writer,"device")
    # Capture contexts are immutable maps; create an independent device instead.
    old=QuantumDevice(:old)
    oldq=add_component!(old,make_qubit(:q;frequency=5.))
    add_recipe!(old,:base,Recipe(old,[oldq]))
    context=QD._encode_record(QD.capture_gate_context(old),writer,"old")
    oldref=only(QD._legacy_fields(context)["recipes"]["entries"])[2]
    # Both contexts deliberately share the old definition identity.
    QD._legacy_fields(definition)["recipes"]["entries"][1][2]=oldref
    records=deepcopy(writer.recipes)
    for record in values(records)
        get(record,"kind",nothing)=="Leaf" && pop!(record["fields"],"input_defaults",nothing)
    end
    historical=Dict("tag"=>"named_tuple","keys"=>["context","recipe"],"values"=>[context,oldref])
    gates=Dict("history"=>historical)
    QD._migrate_recipe_input_records!(definition,records,gates;schema=4)
    reader=QD.RecordReader(records,writer.arrays,Dict(),Set(),Dict(),Dict())
    live=QD._decode_record(definition,reader,"device")
    history=QD._decode_record(gates["history"],reader,"history")
    @test haskey(live.recipes[:base].parameters[string(ν)],9.)
    @test haskey(history.recipe.parameters[string(ν)],5.)
    @test isequal(live.recipes[:base].variables,history.recipe.variables)
end
