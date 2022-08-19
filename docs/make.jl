using EquationsOfStateOfSolids
using Documenter

DocMeta.setdocmeta!(EquationsOfStateOfSolids, :DocTestSetup, :(using EquationsOfStateOfSolids); recursive=true)

makedocs(;
    modules=[EquationsOfStateOfSolids],
    authors="Reno <singularitti@outlook.com>",
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
            "Contributing" => "contributing.md",
        ],
        "API" => [
            "Collections" => "api/collections.md",
            "Finite strains" => "api/finitestrains.md",
            "Fitting" => "api/fitting.md",
            "Portability" => "portability.md",
            "Interoperability" => "interoperability.md",
            "Plotting" => "plotting.md",
        ],
        "Troubleshooting" => "troubleshooting.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/EquationsOfStateOfSolids.jl",
)
