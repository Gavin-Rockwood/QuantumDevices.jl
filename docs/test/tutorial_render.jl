# Isolated end-to-end rendering; no QuantumDevices models are loaded or run.
using Test, Documenter, DocumenterVitepress
include(joinpath(@__DIR__, "..", "tutorial_cache.jl"))

mktempdir() do temporary_root
    root = realpath(temporary_root)
    source = joinpath(root, "fixture.jl")
    write(source, "# # Cached rendering fixture\nvalue = 40\n# Shared state\nvalue + 2\n# Figure\nstruct TinyFigure end\nBase.show(io::IO, ::MIME\"image/svg+xml\", ::TinyFigure) = print(io, \"<svg xmlns='http://www.w3.org/2000/svg' width='2' height='2'><rect width='2' height='2'/></svg>\")\nTinyFigure()\n")
    TutorialCache.tutorial(source, joinpath(root, "src"); repository=root, environment=["fixture"])
    makedocs(; root, sitename="Cache fixture", remotes=nothing, doctest=false,
        format=MarkdownVitepress(; repo="github.com/example/example", devbranch="dev", build_vitepress=false, clean_md_output=false,
            deploy_decision=Documenter.DeployDecision(all_ok=false, repo="example/example", subfolder="dev")),
        pages=["Fixture" => "fixture.md"])
    output = joinpath(root, "build", ".documenter", "fixture.md")
    @test isfile(output)
    markdown = read(output, String)
    @test occursin("42", markdown)
    @test occursin("fixture_assets/page-", markdown)
    @test !occursin("@example", markdown)
    @test any(endswith(".svg"), readdir(joinpath(root, "build", ".documenter", "fixture_assets")))
end
