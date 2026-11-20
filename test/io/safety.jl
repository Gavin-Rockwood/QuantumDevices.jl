@testset "Device bundle failures are explicit and transactional" begin
    QD=QuantumDevices
    root=mktempdir();folder=joinpath(root,"device")
    device=QuantumDevice(:safe)
    save_device(folder,device)
    @test_throws ArgumentError save_device(pwd(),device;overwrite=true)
    @test_throws ArgumentError save_device(dirname(pwd()),device;overwrite=true)
    @test_throws ArgumentError QD._bundle_path(root,"C:\\outside")
    stage=joinpath(root,"stage");mkpath(stage);write(joinpath(stage,"new"),"new")
    original=read(joinpath(folder,"device.yaml"))
    count=Ref(0)
    function failing_rename(from,to)
        count[]+=1
        count[]==2 && error("injected replacement failure")
        mv(from,to)
    end
    @test_throws ErrorException QD._commit_bundle!(stage,folder;overwrite=true,rename=failing_rename)
    @test read(joinpath(folder,"device.yaml"))==original
    @test isdir(stage)
    @test QD._commit_bundle!(stage,folder;overwrite=true)==folder
    @test read(joinpath(folder,"new"),String)=="new"
    # Cyclic records and metadata fail with a useful location rather than recurse forever.
    record=Dict{String,Any}("tag"=>"tuple");record["values"]=[record]
    reader=QD.RecordReader(Dict(),Dict(),Dict(),Set(),Dict(),Dict())
    @test_throws ArgumentError QD._decode_record(record,reader,"cyclic_record")
    array=Any[];push!(array,array)
    @test_throws ArgumentError QD._encode_record(array,QD.RecordWriter(),"cyclic_metadata")
    cyclic_gate=QuantumGate(:cycle,GateRecipe(:model;duration=1.))
    cyclic_gate.metadata=cyclic_gate
    @test_throws ArgumentError QD._encode_record(cyclic_gate,QD.RecordWriter(),"cyclic_gate")
    # Malformed primitive payloads are rejected before Julia object reconstruction.
    save_device(folder,device;overwrite=true)
    manifest=YAML.load_file(joinpath(folder,"device.yaml"))
    payload=joinpath(folder,manifest["payload"])
    QD.JLD2.jldopen(payload,"w") do file;file["bad"]=(;not_primitive=1);end
    manifest["files"][manifest["payload"]]=QD._file_digest(payload)
    YAML.write_file(joinpath(folder,"device.yaml"),manifest)
    @test_throws ArgumentError load_device(folder)
    # Declared path traversal and internal symlinks never resolve outside the bundle.
    save_device(folder,device;overwrite=true)
    manifest=YAML.load_file(joinpath(folder,"device.yaml"));manifest["files"]["../outside"]="invalid"
    YAML.write_file(joinpath(folder,"device.yaml"),manifest)
    @test_throws ArgumentError load_device(folder)
    save_device(folder,device;overwrite=true)
    link=joinpath(folder,"linked");symlink(joinpath(folder,"definitions"),link)
    @test_throws ArgumentError QD._bundle_path(folder,"linked/device.yaml")
end
