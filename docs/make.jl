using EquationsOfStateOfSolids
using Documenter

makedocs(;
    modules=[EquationsOfStateOfSolids],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/blob/{commit}{path}#L{line}",
    sitename="EquationsOfStateOfSolids.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/EquationsOfStateOfSolids.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Manual" => [
            "Collections" => "api/collections.md",
            "Fitting" => "api/fitting.md",
            "Portability" => "portability.md",
            "Interoperability" => "interoperability.md",
            "Plotting" => "plotting.md",
        ],
        "FAQ" => "faq.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/EquationsOfStateOfSolids.jl",
)
