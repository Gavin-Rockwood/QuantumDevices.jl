using QuantumDevices
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(QuantumDevices, :DocTestSetup, :(using QuantumDevices); recursive=true)
const PAGES = [
    "Home" => "index.md",
    "Getting Started" => [
        "Overview" => "getting_started/overview.md",
    ],
    "User Guide" => [
        "Circuits" => [
            "Circuit Elements" => "user_guide/circuits/circuit_elements/circuit_elements.md",
            "Building Circuits" => "user_guide/circuits/building_circuits/building_circuits.md",
        ],
        "Dynamics" => [
            "Floquet Tools" => "user_guide/dynamics/floquet/floquet.md",
        ],
    ],
    "Resources" => [
        "API" => "resources/api.md",
    ],   
]

makedocs(;
    modules = [QuantumDevices],
    authors = "Gavin Rockwood",
    sitename = "QuantumDevices.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/Gavin-Rockwood/QuantumDevices.jl",
        devbranch = "main",
    ),
    remotes = nothing,
    pages = PAGES,
    checkdocs = :none
)

deploydocs(;
    repo = "github.com/Gavin-Rockwood/QuantumDevices.jl",
    target = joinpath(@__DIR__, "build"),
    devbranch = "main",
    branch = "gh-pages",
    push_preview = true,
)
