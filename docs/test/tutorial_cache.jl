using Test
include(joinpath(@__DIR__, "..", "tutorial_cache.jl"))
using .TutorialCache

@testset "Whole tutorial execution cache" begin
    mktempdir() do repository
        source = joinpath(repository, "fixture.jl")
        write(source, "# # Fixture\nx = 20\n# Continued state\nx + 22\n# Output\nprintln(\"fixture output\")\n# Figure\nstruct TinyFigure end\nBase.show(io::IO, ::MIME\"image/svg+xml\", ::TinyFigure) = print(io, \"<svg xmlns='http://www.w3.org/2000/svg' width='2' height='2'><rect width='2' height='2'/></svg>\")\nTinyFigure()\n")
        destination = joinpath(repository, "published")
        environment = ["fixture-dependency:1"]
        count = Ref(0)
        renderer = function (args...)
            count[] += 1
            TutorialCache.render(args...)
        end
        run(;kwargs...) = TutorialCache.tutorial(source, destination; repository,
            environment, renderer, kwargs...)
        @test run() == :executed
        @test count[] == 1
        markdown = read(joinpath(destination, "fixture.md"), String)
        @test occursin("42", markdown)
        @test occursin("fixture output", markdown)
        @test !occursin("@example", markdown)
        @test occursin("fixture_assets/page-", markdown)
        @test any(endswith(".svg"), readdir(joinpath(destination, "fixture_assets")))
        @test !isfile(joinpath(destination, "fixture_assets", "page.md"))
        @test run() == :hit
        @test count[] == 1
        @test run(mode="refresh") == :executed
        @test count[] == 2
        @test run(mode="preview") == :current
        write(source, read(source, String) * "\n# Changed prose\n")
        @test run(mode="preview") == :stale
        @test count[] == 2
        @test occursin("preview: stale", read(joinpath(destination, "fixture.md"), String))
        @test run() == :executed
        before = TutorialCache.fingerprint(source; repository, environment)
        @test before != TutorialCache.fingerprint(source; repository, environment=["fixture-dependency:2"])
        write(joinpath(repository, "data.txt"), "one")
        key1 = TutorialCache.fingerprint(source; repository, environment, inputs=["data.txt"])
        write(joinpath(repository, "data.txt"), "two")
        @test key1 != TutorialCache.fingerprint(source; repository, environment, inputs=["data.txt"])
        @test_throws ErrorException TutorialCache.fingerprint(source; repository, inputs=["missing.txt"])
        # Failure leaves the last good generation usable and is never a fallback.
        failing(args...) = error("fixture failure")
        @test_throws ErrorException TutorialCache.tutorial(source, destination;
            repository, environment, mode="refresh", renderer=failing)
        @test run() == :hit
        @test count[] == 3
        changing = function(args...)
            TutorialCache.render(args...)
            write(source, read(source, String) * "\n# Edited during execution\n")
        end
        @test_throws ErrorException TutorialCache.tutorial(source, destination;
            repository, environment, mode="refresh", renderer=changing)
        @test run(mode="preview") == :stale
        # Restore source to its last successful identity.
        write(source, replace(read(source, String), "\n# Edited during execution\n" => ""))
        @test run() == :hit
        mkpath(joinpath(repository, "src"))
        write(joinpath(repository, "src", "helper.jl"), "value = 1")
        src_key = TutorialCache.fingerprint(source; repository, environment)
        write(joinpath(repository, "src", "helper.jl"), "value = 2")
        @test src_key != TutorialCache.fingerprint(source; repository, environment)
        rm(joinpath(repository, "src"); recursive=true)
        mkpath(joinpath(repository, "docs"))
        write(joinpath(repository, "docs", "Project.toml"), "name = \"Fixture\"\n")
        docs_key = TutorialCache.fingerprint(source; repository, environment)
        write(joinpath(repository, "docs", "Project.toml"), "name = \"ChangedFixture\"\n")
        @test docs_key != TutorialCache.fingerprint(source; repository, environment)
        rm(joinpath(repository, "docs", "Project.toml"))
        # Corrupt/missing assets invalidate an otherwise matching entry.
        namespace = joinpath(repository, "docs", ".execution-cache", "fixture", "full")
        entry = joinpath(namespace, strip(read(joinpath(namespace, "latest"), String)))
        write(joinpath(entry, "output", "page.md"), "corrupted")
        @test run() == :executed
        @test count[] == 4
        emptycache = joinpath(repository, "empty-cache")
        @test TutorialCache.tutorial(source, destination; repository, environment,
            cache=emptycache, mode="preview", renderer=failing) == :missing
        @test occursin("preview: missing", read(joinpath(destination, "fixture.md"), String))
    end
end
