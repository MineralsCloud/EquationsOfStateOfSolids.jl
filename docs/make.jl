using EquationsOfStateOfSolids
using Documenter

DocMeta.setdocmeta!(EquationsOfStateOfSolids, :DocTestSetup, :(using EquationsOfStateOfSolids); recursive=true)

makedocs(;
    modules=[EquationsOfStateOfSolids],
    authors="singularitti <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/blob/{commit}{path}#{line}",
    sitename="EquationsOfStateOfSolids.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/EquationsOfStateOfSolids.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Installation guide" => "installation.md",
            "Plotting" => "plotting.md",
            "Interoperability" => "interoperability.md",
            "Portability" => "portability.md",
        ],
        "API Reference" => [
            "Collections" => "api/collections.md",
            "Finite strains" => "api/finitestrains.md",
            "Fitting" => "api/fitting.md",
        ],
        "Developer Docs" => [
            "Contributing" => "developers/contributing.md",
            "Style Guide" => "developers/style.md",
        ],
        "Troubleshooting" => "troubleshooting.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/EquationsOfStateOfSolids.jl",
)
