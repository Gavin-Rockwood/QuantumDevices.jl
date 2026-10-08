using QuantumDevices
using Documenter
using DocumenterVitepress
using Literate

const DOCS_ROOT = @__DIR__
const BUILD_ROOT = joinpath(DOCS_ROOT, "build")
const MARKDOWN_ONLY = get(ENV, "DOCS_MARKDOWN_ONLY", "false") == "true"
const DEPLOY = get(ENV, "DOCS_DEPLOY", "false") == "true"

# The Julia walkthrough is the source for both the script and tutorial page.
Literate.markdown(
    joinpath(DOCS_ROOT, "src", "tutorials", "transmon_resonator_control.jl"),
    joinpath(DOCS_ROOT, "src", "tutorials");
    name="transmon_resonator_control", execute=false, credit=false,
    codefence="```@example mode3" => "```",
)

Literate.markdown(
    joinpath(DOCS_ROOT, "src", "tutorials", "tunable_coupler_control.jl"),
    joinpath(DOCS_ROOT, "src", "tutorials");
    name="tunable_coupler_control", execute=false, credit=false,
    codefence="```@example coupler" => "```",
    postprocess=content -> rstrip(content) * "\n",
)

function documentation_versions()
    repository = normpath(joinpath(DOCS_ROOT, ".."))
    tags = split(read(`git -C $repository tag --list`, String), '\n'; keepempty=false)
    filter!(tag -> startswith(tag, "v") && tryparse(VersionNumber, tag) !== nothing, tags)
    sort!(tags; by=VersionNumber, rev=true)
    return vcat([tag => tag for tag in tags], ["dev" => "dev"])
end

DocMeta.setdocmeta!(QuantumDevices, :DocTestSetup, :(using QuantumDevices); recursive=true)

# checkdocs catches omitted documented exports; this audit also catches absent docstrings.
undocumented = filter(names(QuantumDevices)) do name
    name != :QuantumDevices && Base.Docs.doc(Base.Docs.Binding(QuantumDevices, name)) === nothing
end
isempty(undocumented) || error("Exported API lacks docstrings: $(join(undocumented, ", "))")

const PAGES = [
    "Home" => "index.md",
    "Getting started" => [
        "Installation and concepts" => "getting_started/overview.md",
        "First model and gate" => "getting_started/quickstart.md",
    ],
    "User guide" => [
        "Symbolic Hamiltonians" => "user_guide/symbolics.md",
        "Components" => "user_guide/components.md",
        "Models and truncation" => "user_guide/models.md",
        "Pulses and flattops" => "user_guide/pulses.md",
        "Gates and evolution" => "user_guide/gates.md",
        "Calibration" => "user_guide/calibration.md",
        "Parameter updates" => "user_guide/parameters.md",
        "State tracking" => "user_guide/tracking.md",
        "Persistence" => "user_guide/persistence.md",
    ],
    "Tutorials" => [
        "Transmon resonator control" => "tutorials/transmon_resonator_control.md",
        "Tunable coupler control" => "tutorials/tunable_coupler_control.md",
    ],
    "Technical explanations" => [
        "Conventions and metrics" => "explanations/conventions.md",
        "Projection and bases" => "explanations/projection.md",
    ],
    "API reference" => [
        "Overview" => "resources/api.md",
        "Symbolics" => "reference/symbolics.md",
        "Components and models" => "reference/models.md",
        "Pulses, gates, calibration" => "reference/gates.md",
        "Tracking and paths" => "reference/tracking.md",
        "Spectral tools" => "reference/spectral_tools.md",
        "Persistence" => "reference/persistence.md",
    ],
    "Development" => [
        "Extensions" => "development/extensions.md",
        "Contributing and builds" => "development/contributing.md",
    ],
]

makedocs(;
    root = DOCS_ROOT,
    modules = [QuantumDevices],
    authors = "Gavin Rockwood",
    sitename = "QuantumDevices.jl",
    format = MarkdownVitepress(
        repo = "github.com/Gavin-Rockwood/QuantumDevices.jl",
        devbranch = "dev",
        description = "Symbolic quantum-device models, shaped controls, and calibrated gates in Julia.",
        build_vitepress = false,
        clean_md_output = false,
        deploy_decision = DEPLOY ? nothing : Documenter.DeployDecision(
            all_ok=false, repo="github.com/Gavin-Rockwood/QuantumDevices.jl", subfolder="dev"),
    ),
    remotes = nothing,
    pages = PAGES,
    checkdocs = :exports,
    doctest = true,
    warnonly = false,
)

if !MARKDOWN_ONLY
    cd(DOCS_ROOT) do
        run(`npm ci --no-audit --no-fund`)
        run(`npm run docs:build`)
    end
    isfile(joinpath(BUILD_ROOT, "final_site", "index.html")) || error("VitePress did not emit index.html")
end

if DEPLOY
    MARKDOWN_ONLY && error("Cannot deploy a Markdown-only build")
    deploydocs(;
        root = DOCS_ROOT,
        repo = "github.com/Gavin-Rockwood/QuantumDevices.jl",
        target = "build/final_site",
        devbranch = "dev",
        branch = "gh-pages",
        push_preview = true,
        versions = documentation_versions(),
    )
end
