import TOML
# Bundle paths never name external files; reject symlinks rather than resolving them implicitly.
function _bundle_path(root,relative,path=relative)
    relative isa AbstractString && !isempty(relative) && !isabspath(relative) && !occursin(r"^[A-Za-z]:",relative) && !occursin('\0',relative) || _ioerror(path,"invalid relative file path")
    parts=split(replace(relative,'\\'=>'/'),'/')
    all(p->!isempty(p) && p!=".." && p!=".",parts) || _ioerror(path,"path traversal is forbidden")
    target=root
    for part in parts
        target=joinpath(target,part)
        islink(target) && _ioerror(path,"symlinks are forbidden")
    end
    return target
end
function _declared_file(root,manifest,relative)
    haskey(manifest["files"],relative) || _ioerror(relative,"file is not declared in manifest")
    return _bundle_path(root,relative)
end
function _write_yaml(filename,record)
    mkpath(dirname(filename))
    YAML.write_file(filename,record)
end
_file_digest(filename)=bytes2hex(SHA.sha256(read(filename)))
function _implementation_digest()
    source=joinpath(pkgdir(QuantumDevices),"src")
    files=sort!([joinpath(root,name) for (root,dirs,names) in walkdir(source) for name in names if endswith(name,".jl") && !endswith(name,"2.jl")])
    return recipe_fingerprint(Dict(relpath(file,source)=>_file_digest(file) for file in files))
end
_environment_fingerprint()=Dict("implementation"=>_implementation_digest(),"julia"=>string(VERSION),"quantumdevices"=>string(Base.pkgversion(QuantumDevices)),
    "manifest"=>isfile(joinpath(pkgdir(QuantumDevices),"Manifest.toml")) ? _file_digest(joinpath(pkgdir(QuantumDevices),"Manifest.toml")) : "none")
function _read_manifest(root;include_native=false)
    islink(root) && _ioerror(root,"bundle root may not be a symlink")
    manifest=YAML.load_file(_bundle_path(root,"device.yaml"))
    get(manifest,"format",nothing)=="quantumdevices.device" && get(manifest,"schema",nothing) in (1,2) || _ioerror(root,"unsupported device format/schema")
    for (relative,digest) in manifest["files"]
        file=_bundle_path(root,relative)
        relative=="checkpoints/native.jld2" && !include_native && continue
        isfile(file) && _file_digest(file)==digest || _ioerror(relative,"missing or corrupt bundle file")
    end
    return manifest
end
_codec_module(type)=type <: Type && !isempty(type.parameters) && first(type.parameters) isa Type ? parentmodule(first(type.parameters)) : parentmodule(type)
function _source_module(type)
    owner=_codec_module(type)
    while parentmodule(owner)!==Main && parentmodule(owner)!==owner
        owner=parentmodule(owner)
    end
    return owner
end
function _module_contains(owner,root)
    while owner!==root && parentmodule(owner)!==owner
        owner=parentmodule(owner)
    end
    return owner===root
end
function _source_descriptor(id,input)
    if input isa AbstractString
        matching=[c for c in values(_IO_CODECS) if c.source_id==id]
        isempty(matching) && _ioerror("code/$id","source has no registered codec")
        module_name=string(nameof(_source_module(first(matching).type)))
        return (;root=dirname(abspath(input)),entrypoint=basename(input),module_name,files=[basename(input)],dependencies=String[])
    end
    all(k->hasproperty(input,k),(:root,:entrypoint,:module_name,:files,:dependencies)) || _ioerror("code/$id","invalid source descriptor")
    return input
end
_source_fingerprint(source)=recipe_fingerprint(Dict("files"=>Dict(k=>bytes2hex(v) for (k,v) in source.bytes),
    "dependencies"=>source.dependencies,"entrypoint"=>source.entrypoint,"module_name"=>source.module_name))
function _check_source_dependencies(dependencies,path)
    for requirement in dependencies
        parts=split(requirement,'=';limit=2);name=String(first(parts))
        Base.isidentifier(name) || _ioerror(path,"invalid dependency name $name")
        location=Base.find_package(name)
        isnothing(location) && _ioerror(path,"required package $name is unavailable; install it explicitly")
        if length(parts)==2
            project=joinpath(dirname(dirname(location)),"Project.toml")
            isfile(project) || _ioerror(path,"cannot verify dependency version for $name")
            actual=get(TOML.parsefile(project),"version",nothing)
            !isnothing(actual) && VersionNumber(actual)==VersionNumber(parts[2]) || _ioerror(path,"dependency $name requires exact version $(parts[2])")
        end
    end
end
function _writer_source!(w,id)
    return get!(w.source_data,id) do
        if haskey(w.source_inputs,id)
            descriptor=_source_descriptor(id,w.source_inputs[id])
            bytes=Dict(String(file)=>read(_bundle_path(abspath(descriptor.root),String(file))) for file in descriptor.files)
            haskey(bytes,String(descriptor.entrypoint)) || _ioerror("code/$id","entrypoint is not declared")
            (;entrypoint=String(descriptor.entrypoint),module_name=String(descriptor.module_name),dependencies=String.(descriptor.dependencies),bytes)
        elseif haskey(_IO_SOURCES,id)
            source=_IO_SOURCES[id]
            (;entrypoint=source.entrypoint,module_name=source.module_name,dependencies=copy(source.dependencies),bytes=Dict(k=>copy(v) for (k,v) in source.bytes))
        else
            _ioerror("code/$id","supply source files for this custom codec")
        end
    end
end
function _capture_sources(w,inputs)
    result=Dict(id=>_writer_source!(w,id) for id in w.sources)
    for (id,source) in result
        _check_source_dependencies(source.dependencies,"code/$id")
        fingerprint=_source_fingerprint(source)
        haskey(_IO_SOURCES,id) && _source_fingerprint(_IO_SOURCES[id])!=fingerprint && _ioerror("code/$id","source ID conflicts with retained provenance")
        module_name=Symbol(source.module_name)
        if haskey(_IO_MODULES,module_name)
            _IO_MODULES[module_name]==fingerprint || _ioerror("code/$id","source conflicts with registered module")
        elseif isdefined(Main,module_name)
            codecs=[c for c in values(_IO_CODECS) if c.source_id==id]
            all(c->_module_contains(_codec_module(c.type),getfield(Main,module_name)),codecs) || _ioerror("code/$id","registered implementations do not belong to supplied module")
            _IO_MODULES[module_name]=fingerprint
        end
    end
    return result
end
function _load_bundle_code!(root,manifest;allow_code)
    for source in manifest["code"]
        id=source["id"];module_name=Symbol(source["module_name"])
        source["entrypoint"] in source["files"] || _ioerror("code/$id","entrypoint must be a declared source file")
        _check_source_dependencies(source["dependencies"],"code/$id")
        bytes=Dict(file=>read(_declared_file(root,manifest,"code/$id/$file")) for file in source["files"])
        fingerprint=_source_fingerprint((;bytes,dependencies=source["dependencies"],entrypoint=source["entrypoint"],module_name=source["module_name"]))
        haskey(_IO_SOURCES,id) && _source_fingerprint(_IO_SOURCES[id])!=fingerprint && _ioerror("code/$id","source ID conflicts with retained provenance")
        if haskey(_IO_MODULES,module_name)
            _IO_MODULES[module_name]==fingerprint || _ioerror("code/$id","module has conflicting source; use a fresh Julia session")

        else
            allow_code || _ioerror("code/$id","bundled Julia source requires allow_code=true")
            isdefined(Main,module_name) && _ioerror("code/$id","Main already contains an unregistered module; use a fresh Julia session")
            Base.include(Main,_declared_file(root,manifest,"code/$id/$(source["entrypoint"])"))
            isdefined(Main,module_name) || _ioerror("code/$id","entrypoint did not define declared module")
            _IO_MODULES[module_name]=fingerprint
        end
        for codec in source["codecs"]
            haskey(_IO_CODECS,codec["id"]) && _IO_CODECS[codec["id"]].version==codec["version"] || _ioerror("code/$id","required codec not registered")
        end
        _IO_SOURCES[id]=(;entrypoint=source["entrypoint"],module_name=String(module_name),dependencies=String.(source["dependencies"]),bytes)
    end
end
