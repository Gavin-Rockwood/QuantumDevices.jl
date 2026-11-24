using YAML,SHA
include("codecs.jl")
include("records.jl")
include("control_maps.jl")
include("symbolics.jl")
include("recipes.jl")
include("generator_migration.jl")
include("parameter_migration.jl")
include("recipe_input_migration.jl")
include("quantum.jl")
include("domain.jl")
include("analysis.jl")
include("manifest.jl")
include("validation.jl")
include("transaction.jl")

# Registry paths and gate names are independent; explicit aliases share one saved gate ID.
function _repository_record(repository,w,gates,identities,path="gates")
    _validate_gate_repository(repository)
    entries=Dict{String,Any}()
    for (name,value) in repository
        if value isa Dict
            entries[String(name)]=Dict("kind"=>"group","entries"=>_repository_record(value,w,gates,identities,"$path/$name"))
        else
            id=get!(identities,value) do
                next="g$(length(gates)+1)"
                gates[next]=_encode_record(value,w,"$path/$name")
                next
            end
            entries[String(name)]=Dict("kind"=>"gate","id"=>id)
        end
    end
    return entries
end
function _restore_repository!(repository,record,gates,path="gates",active=Base.IdSet())
    record in active && _ioerror(path,"cyclic repository manifest")
    push!(active,record)
    try
    for (name,node) in record
        if node["kind"]=="group"
            child=Dict{Symbol,Union{Dict,QuantumGate}}()
            _restore_repository!(child,node["entries"],gates,"$path/$name",active)
            repository[Symbol(name)]=child
        elseif node["kind"]=="gate"
            haskey(gates,node["id"]) || _ioerror(path,"unknown gate ID")
            repository[Symbol(name)]=gates[node["id"]]
        else
            _ioerror(path,"unknown repository node")
        end
    end
    finally
        delete!(active,record)
    end
end
"""
    save_device(folder, device; overwrite=false, code_sources=Dict(), native_checkpoint=false)

Write a complete shareable YAML/JLD2 directory. Recipes are content-interned;
custom callables need registered codecs and explicitly supplied source files.
Only registered gates and attached artifacts are saved. Never builds or evolves.
"""
function save_device(folder::AbstractString,device::QuantumDevice;overwrite=false,code_sources=Dict(),native_checkpoint=false)
    destination=abspath(folder)
    relative_current=relpath(abspath(pwd()),destination)
    (relative_current=="." || !(relative_current==".." || startswith(relative_current,".."*string(Base.Filesystem.path_separator)))) && _ioerror(destination,"destination contains the working directory")
    islink(destination) && _ioerror(destination,"destination symlinks are forbidden")
    ispath(destination) && !overwrite && _ioerror(destination,"destination exists; use overwrite=true")
    ispath(destination) && !isdir(destination) && _ioerror(destination,"destination is not a directory")
    _validate_device_definitions(device)
    w=RecordWriter(;code_sources)
    validate_component_input_defaults(device.components)
    # IO owns encoded records; do not first deepcopy/freeze the whole live device.
    physical=(;name=device.name,components=device.components,interactions=device.interactions,
        recipes=device.recipes,metadata=device.metadata,control_maps=device.control_maps)
    context=_encode_record(physical,w,"device")
    gates=Dict{String,Any}()
    repository=_repository_record(device.gates,w,gates,IdDict())
    sources=_capture_sources(w,code_sources) # Preflight completes before destination changes.
    mkpath(dirname(destination))
    stage=mktempdir(dirname(destination);prefix=".quantumdevices-")
    try
        _write_yaml(joinpath(stage,"definitions/device.yaml"),context)
        for (id,record) in w.recipes
            _write_yaml(joinpath(stage,"recipes/$id.yaml"),record)
        end
        for (id,record) in gates
            _write_yaml(joinpath(stage,"gates/$id.yaml"),record)
        end
        mkpath(joinpath(stage,"payloads"))
        JLD2.jldopen(joinpath(stage,"payloads/arrays.jld2"),"w") do file
            for (id,array) in w.arrays;file[id]=array;end
        end
        source_records=Any[]
        for (id,source) in sources
            for (relative,bytes) in source.bytes
                file=_bundle_path(stage,"code/$id/$relative");mkpath(dirname(file));write(file,bytes)
            end
            codecs=[Dict("id"=>key,"version"=>codec.version) for (key,codec) in _IO_CODECS if codec.source_id==id]
            push!(source_records,Dict("id"=>id,"entrypoint"=>source.entrypoint,"module_name"=>source.module_name,"dependencies"=>source.dependencies,"files"=>sort!(collect(keys(source.bytes))),"codecs"=>codecs))
        end
        if native_checkpoint
            mkpath(joinpath(stage,"checkpoints"))
            JLD2.jldsave(joinpath(stage,"checkpoints/native.jld2");solutions=w.native)
        end
        files=Dict{String,String}()
        for (root,dirs,names) in walkdir(stage), name in names
            file=joinpath(root,name);files[relpath(file,stage)]=_file_digest(file)
        end
        manifest=Dict("format"=>"quantumdevices.device","schema"=>5,"producer"=>_environment_fingerprint(),
            "definition"=>"definitions/device.yaml","recipe_nodes"=>Dict(id=>"recipes/$id.yaml" for id in keys(w.recipes)),
            "gate_nodes"=>Dict(id=>"gates/$id.yaml" for id in keys(gates)),"gates"=>repository,"payload"=>"payloads/arrays.jld2",
            "code"=>source_records,"native_checkpoint"=>native_checkpoint,"files"=>files)
        _write_yaml(joinpath(stage,"device.yaml"),manifest)
        _read_manifest(stage;include_native=native_checkpoint)
        _commit_bundle!(stage,destination;overwrite)
        return destination
    finally
        isdir(stage) && rm(stage;recursive=true)
    end
end
"""
    load_device(folder; allow_code=false, load_native=false, trusted_native=false, migration_report=nothing)

Eagerly restore a complete device without building or simulating. Bundled Julia
source executes in Main only with opt-in. Native checkpoints require independent
trust and exact environment compatibility; portable summaries are the default.
Supply a vector as `migration_report` to collect dropped/unavailable legacy input
notices; otherwise a nonempty migration report gets one informational summary.
"""
function load_device(folder::AbstractString;allow_code=false,load_native=false,trusted_native=false,migration_report=nothing)
    migration_report isa Union{Nothing,AbstractVector} || throw(ArgumentError("migration_report must be a vector or nothing"))
    root=abspath(folder)
    load_native && !trusted_native && _ioerror(root,"native checkpoint requires trusted_native=true")
    manifest=_read_manifest(root;include_native=load_native)
    _load_bundle_code!(root,manifest;allow_code)
    arrays=Dict{String,Any}()
    JLD2.jldopen(_declared_file(root,manifest,manifest["payload"]),"r") do file
        for id in keys(file)
            dataset=JLD2.get_dataset(file,id)
            dataset.datatype == JLD2.h5fieldtype(file,UInt8,UInt8,Val{false}) && isempty(dataset.attributes) ||
                _ioerror("payload/$id","only raw UInt8 datasets are accepted")
            array=JLD2.read_dataset(dataset)
            array isa Vector{UInt8} || _ioerror("payload/$id","nonprimitive numerical payload")
            arrays[id]=array
        end
    end
    native=Dict{String,Any}()
    if load_native
        trusted_native || _ioerror(root,"native checkpoint requires trusted_native=true")
        manifest["native_checkpoint"] || _ioerror(root,"bundle has no native checkpoint")
        manifest["producer"]==_environment_fingerprint() || _ioerror(root,"native checkpoint environment is incompatible")
        native=JLD2.load(_declared_file(root,manifest,"checkpoints/native.jld2"),"solutions")
    end
    recipes=Dict(id=>YAML.load_file(_declared_file(root,manifest,file)) for (id,file) in manifest["recipe_nodes"])
    all(recipe_fingerprint(record)==id for (id,record) in recipes) || _ioerror(root,"recipe content ID mismatch")
    definition=YAML.load_file(_declared_file(root,manifest,manifest["definition"]))
    gate_records=Dict(id=>YAML.load_file(_declared_file(root,manifest,file)) for (id,file) in manifest["gate_nodes"])
    if manifest["schema"]==1
        foreach(_migrate_control_map_records!,values(recipes))
        _migrate_control_map_records!(definition)
        foreach(_migrate_control_map_records!,values(gate_records))
    end
    _migrate_operator_records!(definition,recipes,gate_records;schema=manifest["schema"])
    report = NamedTuple[]
    _migrate_parameter_records!(definition,recipes,gate_records;schema=manifest["schema"],report)
    _migrate_recipe_input_records!(definition,recipes,gate_records;schema=manifest["schema"])
    isnothing(migration_report) || append!(migration_report,report)
    isnothing(migration_report) && !isempty(report) && @info "Device parameter migration" notices=length(report)
    r=RecordReader(recipes,arrays,Dict(),Set(),Dict(),native)
    context=Base.invokelatest(_decode_record,definition,r,"device")
    manifest["schema"] in (2,3,4,5) && !hasproperty(context,:control_maps) && _ioerror("device","schema 2/3/4 requires the control_maps registry field")
    device=Base.invokelatest(QuantumDevice,context)
    gates=Dict(id=>Base.invokelatest(_decode_record,record,r,"gates/$id") for (id,record) in gate_records)
    for (id,gate) in gates;Base.invokelatest(_validate_stored_gate,gate,"gates/$id");end
    _restore_repository!(device.gates,manifest["gates"],gates)
    _validate_device_definitions(device)
    return device
end
