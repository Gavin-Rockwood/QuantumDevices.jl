using YAML, SHA

@testset "current recipe values and legacy default materialization" begin
    QD=QuantumDevices
    d=QuantumDevice(:legacy_inputs)
    q=add_component!(d,make_qubit(:q;frequency=5.))
    g=Symbolics.variable(:legacy_g)
    retired=Symbolics.variable(:retired_default)
    recipe=Recipe(d,[q];interactions=[Interaction(:shift,g*q.operators["z"])])
    add_recipe!(d,:model,recipe)
    folder=joinpath(mktempdir(),"legacy")
    save_device(folder,d)
    manifest=YAML.load_file(joinpath(folder,"device.yaml"))
    @test manifest["schema"]==5
    definition=YAML.load_file(joinpath(folder,manifest["definition"]))
    fields=QD._legacy_fields(definition)
    writer=QD.RecordWriter()
    fields["explicit_parameters"]=QD._encode_record(Dict(g=>.12,retired=>3.),writer,"legacy")
    QD._set_legacy_fields!(definition,fields)
    YAML.write_file(joinpath(folder,manifest["definition"]),definition)
    manifest["schema"]=3
    manifest["files"][manifest["definition"]]=bytes2hex(SHA.sha256(read(joinpath(folder,manifest["definition"]))))
    YAML.write_file(joinpath(folder,"device.yaml"),manifest)
    original_bytes=read(joinpath(folder,"device.yaml"))
    report=NamedTuple[]
    restored=load_device(folder;migration_report=report)
    @test !hasproperty(restored,:explicit_parameters)
    @test restored.recipes[:model].overrides[g]==.12
    @test build(restored,:model).resolved_parameters[g]==.12
    @test haskey(restored.parameters[string(g)],.12)
    @test !haskey(restored.variables,string(retired))
    @test any(item->item.variable==string(retired),report)
    @test read(joinpath(folder,"device.yaml"))==original_bytes
    save_device(joinpath(dirname(folder),"current"),restored)
    current=load_device(joinpath(dirname(folder),"current"))
    @test current.parameters==restored.parameters
end

@testset "historical parameter migration authority and gate-owned copies" begin
    QD=QuantumDevices
    d=QuantumDevice(:contexts)
    q=add_component!(d,make_qubit(:q;frequency=5.))
    ν=q.variables["ν"]
    A=Symbolics.variable(:migration_drive)
    base=Recipe(d,[q]);add_recipe!(d,:base,base)
    endpoint=draft(base);endpoint.overrides[A]=.1
    term=Interaction(:drive,A*q.operators["x"])
    gate=GateRecipe(endpoint;terms=[term],controls=[Control(:drive,A,Pulse(Constant(;duration=.01);amplitude=.1))])
    prepared=build(d,gate;parameters=Dict(A=>.4,ν=>5.4))
    problem=CalibrationProblem(prepared,["drive/amplitude"],Matrix{ComplexF64}(I,2,2);
        states=[(0,),(1,)],lb=[0.],ub=[1.])
    result=CalibrationResult(problem,(;u=[.2],retcode=QD.SciMLBase.ReturnCode.MaxIters))
    q.parameters["ν"]=9. # The live context is intentionally different from the snapshot.
    writer=QD.RecordWriter()
    definition=QD._encode_record(QD.capture_gate_context(d),writer,"device")
    historical=QD._encode_record(prepared,writer,"prepared")
    named=QD._encode_record(GateRecipe(:base;terms=[term],duration=.01),writer,"named")
    model=QD._encode_record(prepared.model,writer,"model")
    calibration=QD._encode_record(result,writer,"calibration")
    # Old device defaults and captured defaults are separate owning contexts.
    fields=QD._legacy_fields(definition)
    inactive=Symbolics.variable(:inactive_old_override)
    fields["explicit_parameters"]=QD._encode_record(Dict(A=>.7,inactive=>3.),writer,"defaults")
    base_ref=only(fields["recipes"]["entries"])[2]
    writer.recipes[base_ref["id"]]["fields"]["parameters"]=QD._encode_record(Dict(inactive=>4.),writer,"inactive")
    QD._set_legacy_fields!(definition,fields)
    context=historical["fields"]["context"]
    captured=QD._legacy_fields(context)
    captured["explicit_parameters"]=QD._encode_record(Dict(A=>.2,inactive=>3.),writer,"captured")
    QD._set_legacy_fields!(context,captured)
    calibration_prepared=QD._legacy_fields(calibration["fields"]["provenance"])["prepared"]
    calibration_context=calibration_prepared["fields"]["context"]
    calibration_fields=QD._legacy_fields(calibration_context)
    calibration_fields["explicit_parameters"]=QD._encode_record(Dict(A=>.2,inactive=>3.),writer,"captured calibration")
    QD._set_legacy_fields!(calibration_context,calibration_fields)
    # Remove the old explicit assignment so the removed default needs materializing.
    endpoint_ref=historical["fields"]["recipe"]["fields"]["model_recipe"]
    record=writer.recipes[endpoint_ref["id"]]
    record["fields"]["parameters"]=QD._encode_record(Dict{Symbolics.Num,Number}(),writer,"empty")
    records=deepcopy(writer.recipes)
    gates=Dict("historical"=>historical,"named"=>named,"model"=>model,"calibration"=>calibration)
    report=NamedTuple[]
    QD._migrate_parameter_records!(definition,records,gates;schema=3,report)
    reader=QD.RecordReader(records,writer.arrays,Dict(),Set(),Dict(),Dict())
    restored=QD._decode_record(gates["historical"],reader,"historical")
    @test restored.parameters[A]==.4
    @test restored.parameters[ν]==5.4
    @test restored.recipe.model_recipe.overrides[A]==.4
    @test restored.model.recipe.resolved_parameters[ν]==5.4
    @test !hasproperty(restored.context,:explicit_parameters)
    named_gate=QD._decode_record(gates["named"],reader,"named")
    @test named_gate.model_recipe isa Recipe
    @test named_gate.model_recipe.overrides[A]==.7
    live=QD._decode_record(definition,reader,"device")
    @test !haskey(live.recipes[:base].overrides,A)
    @test !haskey(live.recipes[:base].overrides,inactive)
    @test any(item->get(item,:reason,nothing)===:inactive_override,report)
    @test live.components[:q].parameters["ν"]==9.
    old_model=QD._decode_record(gates["model"],reader,"model")
    @test old_model.recipe.resolved_parameters[ν]==5.4
    @test !isequal(old_model.recipe.resolved_parameters[ν],live.components[:q].parameters["ν"])
    calibrated=QD._decode_record(gates["calibration"],reader,"calibration")
    @test calibrated.recipe.model_recipe.overrides[A]==.4
    @test calibrated.provenance.prepared.parameters[A]==.4
    @test calibrated.provenance.prepared.model.recipe.resolved_parameters[ν]==5.4
    @test CalibrationProblem(calibrated).initial_values==(.2,)
    @test make_hamiltonian(calibrated.provenance.prepared) isa QT.QobjEvo
end

@testset "legacy interaction defaults become owning recipe inputs" begin
    QD=QuantumDevices
    device=QuantumDevice(:interaction_defaults)
    q=add_component!(device,make_qubit(:q))
    coupling=Symbolics.variable(:legacy_interaction_input)
    term=add_interaction!(device,Interaction(:shift,coupling*q.operators["z"]))
    add_recipe!(device,:model,Recipe(device,[q];interactions=[term]))
    writer=QD.RecordWriter()
    definition=QD._encode_record(QD.capture_gate_context(device),writer,"device")
    records=deepcopy(writer.recipes)
    defaults=QD._encode_record(Dict(string(coupling)=>.13),writer,"legacy defaults")
    # Recreate both stores in a schema-3 registry and its captured recipe-local term.
    fields=QD._legacy_fields(definition)
    only(fields["interactions"]["entries"])[2]["fields"]["parameters"]=deepcopy(defaults)
    for record in values(records)
        record["kind"]=="Leaf" && continue
        for interaction in record["fields"]["interactions"]["values"]
            interaction["fields"]["parameters"]=deepcopy(defaults)
        end
    end
    report=NamedTuple[]
    QD._migrate_parameter_records!(definition,records,Dict{String,Any}();schema=3,report)
    reader=QD.RecordReader(records,writer.arrays,Dict(),Set(),Dict(),Dict())
    restored=QuantumDevice(QD._decode_record(definition,reader,"device"))
    @test restored.recipes[:model].overrides[coupling]==.13
    @test build(restored,:model).resolved_parameters[coupling]==.13
    @test haskey(restored.parameters[string(coupling)],.13)
    @test !hasproperty(restored.interactions[:shift],:parameters)
    @test !hasproperty(only(restored.recipes[:model].interactions),:parameters)
    @test isempty(report)
end
