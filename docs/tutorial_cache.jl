module TutorialCache
using SHA, TOML, Pkg
using Literate

const FORMAT = "whole-tutorial-v1"
full_source(source) = replace(source, r"(?ms)^# cheap-docs-begin\n.*?^# cheap-docs-end\n" => "")
cheap_source(source) = replace(source, r"(?m)^# cheap-docs-(?:begin|end)\n" => "")

"""Content fingerprint, using repository-relative names and resolved dependency identities.
Declare every extra included source/data file in `inputs` (relative to `repository`)."""
function fingerprint(source; repository, inputs=String[], settings="full", environment=nothing)
    files = [relpath(source, repository); inputs]
    for directory in ("src",)
        root = joinpath(repository, directory)
        isdir(root) || continue
        for (dir, _, names) in walkdir(root), name in names
            push!(files, relpath(joinpath(dir, name), repository))
        end
    end
    append!(files, ["Project.toml", "Manifest.toml", "docs/Project.toml", "docs/Manifest.toml"])
    io = IOBuffer()
    print(io, FORMAT, '\n', VERSION, '\n', settings, '\n')
    print(io, Sys.ARCH, '\n', Threads.nthreads(), '\n')
    for key in ("JULIA_NUM_THREADS", "JULIA_CPU_TARGET", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")
        print(io, key, '=', get(ENV, key, "<unset>"), '\n')
    end
    for file in sort!(unique(files))
        print(io, file, '\0')
        path = joinpath(repository, file)
        if isfile(path)
            if endswith(file, "Project.toml") || endswith(file, "Manifest.toml")
                # Local path dependencies inside this checkout remain portable.
                print(io, replace(read(path, String), repository => "<repository>"))
            else
                write(io, read(path))
            end
        elseif file in inputs || file == relpath(source, repository)
            error("Missing declared tutorial input: $file")
        else
            print(io, "<absent>")
        end
        print(io, '\0')
    end
    # Include helper bytes and installed dependency versions/tree identities. Avoid
    # absolute checkout/depot paths so identical environments can share a cache.
    write(io, read(@__FILE__))
    dependencies = environment === nothing ? sort!([
        string(uuid, ':', dep.name, ':', dep.version, ':', dep.tree_hash)
        for (uuid, dep) in Pkg.dependencies()]) : environment
    print(io, repr(dependencies))
    return bytes2hex(sha256(take!(io)))
end

function inventory(directory)
    result = Dict{String,String}()
    for (root, _, files) in walkdir(directory), file in files
        path = joinpath(root, file)
        result[relpath(path, directory)] = bytes2hex(sha256(read(path)))
    end
    result
end

function valid_entry(entry, key=nothing)
    try
        metadata = TOML.parsefile(joinpath(entry, "metadata.toml"))
        metadata["format"] == FORMAT || return false
        key === nothing || metadata["fingerprint"] == key || return false
        files = inventory(joinpath(entry, "output"))
        !isempty(files) && files == metadata["files"] && isfile(joinpath(entry, "output", "page.md"))
    catch
        false
    end
end

function render(source, output, name)
    Literate.markdown(source, output; name="page", execute=true,
        flavor=Literate.DocumenterFlavor(), credit=false,
        preprocess=full_source, codefence="```julia" => "```",
        continue_on_error=false)
end

function publish(entry, destination, name; banner="")
    output = joinpath(entry, "output")
    assets = joinpath(destination, name * "_assets")
    mkpath(destination)
    # Own only this tutorial's namespace; Documenter copies its contents normally.
    isdir(assets) && rm(assets; recursive=true)
    mkpath(assets)
    for file in readdir(output)
        file == "page.md" && continue
        cp(joinpath(output, file), joinpath(assets, file))
    end
    markdown = read(joinpath(output, "page.md"), String)
    # Literate emits local image links. Preserve their names and relative paths.
    markdown = replace(markdown, r"!\[([^\]]*)\]\(([^)]+)\)" => matchtext -> begin
        matchvalue = match(r"!\[([^\]]*)\]\(([^)]+)\)", matchtext)
        target = matchvalue.captures[2]
        isfile(joinpath(output, target)) ? "![$(matchvalue.captures[1])]($(name)_assets/$target)" : matchtext
    end)
    write(joinpath(destination, name * ".md"), banner * markdown)
end

function resolve_entry(namespace, pointer)
    isfile(pointer) || return ""
    generation = strip(read(pointer, String))
    occursin(r"^[0-9a-f]{64}--pending-[A-Za-z0-9]+$", generation) || return ""
    return joinpath(namespace, generation)
end

function write_pointer(path, generation)
    mkpath(dirname(path))
    temporary, io = mktemp(dirname(path))
    try
        write(io, generation)
        close(io)
        Base.Filesystem.rename(temporary, path)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
    end
end

"""Generate/reuse a complete tutorial. Preview never executes and returns its status."""
function tutorial(source, destination; name=splitext(basename(source))[1], repository,
                  cache=joinpath(repository, "docs", ".execution-cache"),
                  mode="changed", inputs=String[], environment=nothing, renderer=render)
    mode in ("changed", "refresh", "preview") || error("Invalid DOCS_EXECUTION=$mode")
    key = fingerprint(source; repository, inputs, environment)
    namespace = joinpath(cache, name, "full")
    index = joinpath(namespace, "keys", key)
    entry = resolve_entry(namespace, index)
    hit = valid_entry(entry, key)
    if mode == "preview"
        status = hit ? "current" : "missing"
        if !hit && isfile(joinpath(namespace, "latest"))
            previous = resolve_entry(namespace, joinpath(namespace, "latest"))
            if valid_entry(previous)
                entry = previous
                status = "stale"
            end
        end
        banner = "> **Tutorial execution preview: $status.** No examples or doctests run in this build; deployment is disabled.\n\n"
        if status == "missing"
            mkpath(destination)
            Literate.markdown(source, destination; name, execute=false, documenter=false,
                credit=false, preprocess=full_source, codefence="```julia" => "```",
                postprocess=text -> banner * text)
        else
            publish(entry, destination, name; banner)
        end
        @info "Tutorial preview" name status
        return Symbol(status)
    end
    if hit && mode == "changed"
        @info "Tutorial cache hit" name
    else
        @info "Tutorial execution starting" name reason=(mode == "refresh" ? "refresh" : "missing/stale")
        mkpath(namespace)
        temporary = mktempdir(namespace; prefix="pending-")
        try
            output = joinpath(temporary, "output")
            mkpath(output)
            renderer(source, output, name)
            key == fingerprint(source; repository, inputs, environment) ||
                error("Tutorial inputs changed during execution: $name; result was not cached")
            isfile(joinpath(output, "page.md")) || error("Tutorial renderer did not produce page.md")
            open(joinpath(temporary, "metadata.toml"), "w") do io
                TOML.print(io, Dict("format" => FORMAT, "fingerprint" => key, "files" => inventory(output)))
            end
            # Immutable generations plus atomically replaced small pointers keep
            # the previous successful cache reachable even during refresh/failure.
            generation = key * "--" * basename(temporary)
            entry = joinpath(namespace, generation)
            mv(temporary, entry)
            write_pointer(index, generation)
            write_pointer(joinpath(namespace, "latest"), generation)
            @info "Tutorial execution complete" name
        finally
            isdir(temporary) && rm(temporary; recursive=true)
        end
    end
    valid_entry(entry, key) || error("Invalid tutorial cache: $name")
    publish(entry, destination, name)
    return hit && mode == "changed" ? :hit : :executed
end
end
