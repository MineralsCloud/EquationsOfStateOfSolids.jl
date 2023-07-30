using EquationsOfStateOfSolids
using Documenter

DocMeta.setdocmeta!(
    EquationsOfStateOfSolids,
    :DocTestSetup,
    quote
        using EquationsOfStateOfSolids#=, EquationsOfStateOfSolids.Fitting=#
        using EquationsOfStateOfSolids.FiniteStrains:
            EulerianStrainFromVolume, VolumeFromEulerianStrain
        using Unitful, UnitfulAtomic
    end;
    recursive=true,
)

makedocs(;
    modules=[EquationsOfStateOfSolids],
    authors="singularitti <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/blob/{commit}{path}#{line}",
    sitename="EquationsOfStateOfSolids.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/EquationsOfStateOfSolids.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Installation Guide" => "man/installation.md",
            "Definitions and conventions" => "man/definitions.md",
            "Examples" => "man/examples.md",
            "Troubleshooting" => "man/troubleshooting.md",
        ],
        "Reference" => Any[
            "Public API" => map(
                s -> "lib/public/$(s)",
                sort(readdir(joinpath(@__DIR__, "src/lib/public")))
            ),
            "Internals" => map(
                s -> "lib/internals/$(s)",
                sort(readdir(joinpath(@__DIR__, "src/lib/internals")))
            ),
        ],
        "Developer Docs" => [
            "Contributing" => "developers/contributing.md",
            "Style Guide" => "developers/style-guide.md",
            "Design Principles" => "developers/design-principles.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/EquationsOfStateOfSolids.jl",
    devbranch="main",
)
