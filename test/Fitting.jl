module Fitting

using Test
using Unitful, UnitfulAtomic
using YAML

using EquationsOfStateOfSolids.Collections:
    Parameters, BirchMurnaghan3rd, Murnaghan, PoirierTarantola3rd, Vinet, EnergyEOS
using EquationsOfStateOfSolids.Fitting: nonlinfit, linfit

# Do not export! Only for internal use!
_getfields(x) = (getfield(x, i) for i in 1:nfields(x))
_isapprox(a::T, b::T; kwargs...) where {T<:Parameters} =
    isapprox(ustrip.(_getfields(a)), ustrip.(_getfields(b)); kwargs...)

# Data from https://github.com/materialsproject/pymatgen/blob/19c4d98/pymatgen/analysis/tests/test_eos.py#L17-L73
@testset "Test data from Pymatgen" begin
    data = open("test/data/si.yml", "r") do io
        YAML.load(io)
    end
    volumes, energies = data["volume"], data["energy"]
    @test _isapprox(
        nonlinfit(EnergyEOS(BirchMurnaghan3rd(40, 0.5, 4, 0)), volumes, energies),
        BirchMurnaghan3rd(
            40.98926572528106,
            0.5369258245417454,
            4.178644235500821,
            -10.842803908240892,
        ),
    )
    @test _isapprox(
        linfit(EnergyEOS(BirchMurnaghan3rd(40, 0.5, 4, 0)), volumes, energies),
        BirchMurnaghan3rd(
            40.98926572528106,
            0.5369258245417454,
            4.178644235500821,
            -10.842803908240892,
        ),
    )
    @test _isapprox(
        nonlinfit(EnergyEOS(Murnaghan(41, 0.5, 4, 0)), volumes, energies),
        Murnaghan(
            41.13757930387086,
            0.5144967693786603,
            3.9123862262572264,
            -10.836794514626673,
        ),
    )
    @test _isapprox(
        nonlinfit(EnergyEOS(PoirierTarantola3rd(41, 0.5, 4, 0)), volumes, energies),
        PoirierTarantola3rd(
            40.86770643373908,
            0.5667729960804602,
            4.331688936974368,
            -10.851486685041658,
        ),
    )
    @test _isapprox(
        linfit(EnergyEOS(PoirierTarantola3rd(41, 0.5, 4, 0)), volumes, energies),
        PoirierTarantola3rd(
            40.86770643373908,
            0.5667729960804602,
            4.331688936974368,
            -10.851486685041658,
        ),
    )
    @test _isapprox(
        nonlinfit(EnergyEOS(Vinet(41, 0.5, 4, 0)), volumes, energies),
        Vinet(
            40.916875663779784,
            0.5493839425156859,
            4.3051929654936885,
            -10.846160810560756,
        ),
    )
end

# Data from https://github.com/materialsproject/pymatgen/blob/19c4d98/pymatgen/analysis/tests/test_eos.py#L92-L167
@testset "Test Mg dataset" begin
    data = open("test/data/mp153.yml", "r") do io
        YAML.load(io)
    end
    volumes, energies, known_energies_vinet =
        data["volume"], data["energy"], data["known_energy_vinet"]
    fitted_eos = nonlinfit(EnergyEOS(Vinet(23, 0.5, 4, -2)), volumes, energies)
    @test _isapprox(
        fitted_eos,
        Vinet(22.957645593, 0.22570911414, 4.06054339, -1.59442926062),
    )
    @test isapprox(map(EnergyEOS(fitted_eos), volumes), known_energies_vinet; atol = 1e-6)
end
    )
    @test isapprox(map(EnergyEOS(fitted_eos), volumes), known_energies_vinet; atol = 1e-5)
end

end
