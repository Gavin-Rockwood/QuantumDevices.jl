using QuantumDevices
using Documenter
using DocumenterVitepress
using Literate

const DOCS_ROOT = @__DIR__
const BUILD_ROOT = joinpath(DOCS_ROOT, "build")
const MARKDOWN_ONLY = get(ENV, "DOCS_MARKDOWN_ONLY", "false") == "true"
const DEPLOY = get(ENV, "DOCS_DEPLOY", "false") == "true"
const EXECUTION = get(ENV, "DOCS_EXECUTION", "changed")
EXECUTION in ("changed", "refresh", "preview") || error("DOCS_EXECUTION must be changed, refresh, or preview")
const PREVIEW = EXECUTION == "preview"
const DRAFT = PREVIEW || get(ENV, "DOCS_DRAFT", "false") == "true"
const RUN_TUTORIALS = get(ENV, "DOCS_RUN_TUTORIALS", "false") == "true"
EXECUTION == "refresh" && !RUN_TUTORIALS && error("DOCS_EXECUTION=refresh requires DOCS_RUN_TUTORIALS=true")
DRAFT && RUN_TUTORIALS && !PREVIEW && error("DOCS_RUN_TUTORIALS=true requires DOCS_DRAFT=false; draft builds do not execute examples")
DRAFT && DEPLOY && error("Draft documentation cannot be deployed")

include("tutorial_cache.jl")
# Extra includes/data must be declared here, relative to the repository root.
const TUTORIAL_INPUTS = Dict{String,Vector{String}}()
function generate_tutorials()
    sources = joinpath(DOCS_ROOT, "examples", "tutorials")
    destination = joinpath(DOCS_ROOT, "src", "tutorials")
    for name in ("transmon_resonator_control", "tunable_coupler_control")
        source = joinpath(sources, name * ".jl")
        if RUN_TUTORIALS || PREVIEW
            TutorialCache.tutorial(source, destination; name,
                repository=dirname(DOCS_ROOT), mode=EXECUTION,
                inputs=get(TUTORIAL_INPUTS, name, String[]))
        else
            Literate.markdown(source, destination; name, execute=false,
                documenter=false, credit=false, preprocess=TutorialCache.cheap_source,
                codefence="```julia" => "```")
        end
    end
end
generate_tutorials()

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
    name != :QuantumDevices && !Base.Docs.hasdoc(QuantumDevices, name)
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
        "Gate recipes" => "user_guide/gate_recipes.md",
        "Gates and evolution" => "user_guide/gates.md",
        "Calibration" => "user_guide/calibration.md",
        "Floquet and Stark shifts" => "user_guide/floquet.md",
        "Parameter updates" => "user_guide/parameters.md",
        "State tracking" => "user_guide/tracking.md",
        "Interaction-ramp labels" => "user_guide/state_labeling.md",
        "Persistence" => "user_guide/persistence.md",
    ],
    "Tutorials" => [
        "Transmon–resonator control" => "tutorials/transmon_resonator_control.md",
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
    pagesonly = DRAFT,
    clean = !DRAFT,
    checkdocs = :exports,
    draft = DRAFT,
    doctest = !DRAFT,
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
