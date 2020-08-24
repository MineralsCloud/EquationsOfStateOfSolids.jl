module Fitting

using Test
using Unitful, UnitfulAtomic
using YAML

using EquationsOfStateOfSolids.Collections:
    Parameters,
    Murnaghan,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    BirchMurnaghan4th,
    PoirierTarantola3rd,
    Vinet,
    EnergyEOS
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
        nonlinfit(EnergyEOS(BirchMurnaghan3rd(40, 0.5, 4)), volumes, energies),
        BirchMurnaghan3rd(
            40.98926572528106,
            0.5369258245417454,
            4.178644235500821,
            -10.842803908240892,
        ),
    )
    @test _isapprox(
        linfit(EnergyEOS(BirchMurnaghan3rd(40, 0.5, 4)), volumes, energies),
        BirchMurnaghan3rd(
            40.98926572528106,
            0.5369258245417454,
            4.178644235500821,
            -10.842803908240892,
        ),
    )
    @test _isapprox(
        nonlinfit(EnergyEOS(Murnaghan(41, 0.5, 4)), volumes, energies),
        Murnaghan(
            41.13757930387086,
            0.5144967693786603,
            3.9123862262572264,
            -10.836794514626673,
        ),
    )
    @test _isapprox(
        nonlinfit(EnergyEOS(PoirierTarantola3rd(41, 0.5, 4)), volumes, energies),
        PoirierTarantola3rd(
            40.86770643373908,
            0.5667729960804602,
            4.331688936974368,
            -10.851486685041658,
        ),
    )
    @test _isapprox(
        linfit(EnergyEOS(PoirierTarantola3rd(41, 0.5, 4)), volumes, energies),
        PoirierTarantola3rd(
            40.86770643373908,
            0.5667729960804602,
            4.331688936974368,
            -10.851486685041658,
        ),
    )
    @test _isapprox(
        nonlinfit(EnergyEOS(Vinet(41, 0.5, 4)), volumes, energies),
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
    @testset "without unit" begin
        volumes, energies, known_energies_vinet =
            data["volume"], data["energy"], data["known_energy_vinet"]
        fitted_eos = nonlinfit(EnergyEOS(Vinet(23, 0.5, 4)), volumes, energies)
        @test _isapprox(
            fitted_eos,
            Vinet(22.957645593, 0.22570911414, 4.06054339, -1.59442926062),
        )
        @test isapprox(
            map(EnergyEOS(fitted_eos), volumes),
            known_energies_vinet;
            atol = 1e-6,
        )
    end

    @testset "with units" begin
        volumes, energies, known_energies_vinet = data["volume"] * u"angstrom^3",
        data["energy"] * u"eV",
        data["known_energy_vinet"] * u"eV"
        fitted_eos =
            nonlinfit(EnergyEOS(Vinet(23u"angstrom^3", 36.16u"GPa", 4)), volumes, energies)
        @test _isapprox(
            fitted_eos,
            Vinet(
                22.957645593u"angstrom^3",
                36.16258657649159u"GPa",  # https://github.com/materialsproject/pymatgen/blob/19c4d98/pymatgen/analysis/tests/test_eos.py#L181
                4.06054339,
                -1.59442926062u"eV",
            ),
        )
        @test isapprox(
            map(EnergyEOS(fitted_eos), volumes),
            known_energies_vinet;
            atol = 1e-6u"eV",
        )
    end
end

# Data from https://github.com/materialsproject/pymatgen/blob/19c4d98/pymatgen/analysis/tests/test_eos.py#L185-L260
@testset "Test Si dataset" begin
    data = open("test/data/mp149.yml", "r") do io
        YAML.load(io)
    end
    @testset "without unit" begin
        volumes, energies, known_energies_vinet =
            data["volume"], data["energy"], data["known_energy_vinet"]
        fitted_eos = nonlinfit(EnergyEOS(Vinet(20, 0.5, 4)), volumes, energies)
        @test _isapprox(
            fitted_eos,
            Vinet(20.446696754, 0.55166385214, 4.32437391, -5.42496338987),
        )
        @test isapprox(
            map(EnergyEOS(fitted_eos), volumes),
            known_energies_vinet;
            atol = 1e-5,
        )
    end

    @testset "with units" begin
        volumes, energies, known_energies_vinet = data["volume"] * u"angstrom^3",
        data["energy"] * u"eV",
        data["known_energy_vinet"] * u"eV"
        fitted_eos =
            nonlinfit(EnergyEOS(Vinet(20u"angstrom^3", 88.39u"GPa", 4)), volumes, energies)
        @test _isapprox(
            fitted_eos,
            Vinet(
                20.446696754u"angstrom^3",
                88.38629264585195u"GPa",  # https://github.com/materialsproject/pymatgen/blob/19c4d98/pymatgen/analysis/tests/test_eos.py#L274
                4.32437391,
                -5.42496338987u"eV",
            ),
        )
        @test isapprox(
            map(EnergyEOS(fitted_eos), volumes),
            known_energies_vinet;
            atol = 1e-6u"eV",
        )
    end
end

# Data from https://github.com/materialsproject/pymatgen/blob/19c4d98/pymatgen/analysis/tests/test_eos.py#L278-L353
@testset "Test Ti dataset" begin
    data = open("test/data/mp72.yml", "r") do io
        YAML.load(io)
    end
    @testset "without unit" begin
        volumes, energies, known_energies_vinet =
            data["volume"], data["energy"], data["known_energy_vinet"]
        fitted_eos = nonlinfit(EnergyEOS(Vinet(17, 0.5, 4)), volumes, energies)
        @test _isapprox(
            fitted_eos,
            Vinet(17.1322302613, 0.70297662247, 3.638807756, -7.89741495912),
        )
        @test isapprox(
            map(EnergyEOS(fitted_eos), volumes),
            known_energies_vinet;
            atol = 1e-5,
        )
    end

    @testset "with units" begin
        volumes, energies, known_energies_vinet = data["volume"] * u"angstrom^3",
        data["energy"] * u"eV",
        data["known_energy_vinet"] * u"eV"
        fitted_eos =
            nonlinfit(EnergyEOS(Vinet(17u"angstrom^3", 112.63u"GPa", 4)), volumes, energies)
        @test _isapprox(
            fitted_eos,
            Vinet(
                17.1322302613u"angstrom^3",
                112.62927094503254u"GPa",  # https://github.com/materialsproject/pymatgen/blob/19c4d98/pymatgen/analysis/tests/test_eos.py#L367
                3.638807756,
                -7.89741495912u"eV",
            ),
        )
        @test isapprox(
            map(EnergyEOS(fitted_eos), volumes),
            known_energies_vinet;
            atol = 1e-7u"eV",
        )
    end
end

@testset "Test `w2k-lda-na.dat` from `Gibbs2`" begin
    data = open("test/data/w2k-lda-na.yml", "r") do io
        YAML.load(io)
    end
    @testset "without unit" begin
        volumes, energies = data["volume"], data["energy"]
        # See https://github.com/aoterodelaroza/asturfit/blob/4de9b41/test/test03.out#L117-L122
        @test _isapprox(
            nonlinfit(EnergyEOS(Murnaghan(224, 6e-4, 4)), volumes, energies),
            Murnaghan(224.501825, 6.047952315268776e-4, 3.723835, -323.417686);
            atol = 1e-5,
        )
        # No reference data, I run on my computer.
        @test _isapprox(
            nonlinfit(EnergyEOS(BirchMurnaghan2nd(224, 6e-4)), volumes, energies),
            BirchMurnaghan2nd(223.7192539523166, 6.268341030294978e-4, -323.4177121144877);
            atol = 1e-8,
        )
        # See https://github.com/aoterodelaroza/asturfit/blob/4de9b41/test/test03.out#L15-L20
        @test _isapprox(
            nonlinfit(EnergyEOS(BirchMurnaghan3rd(224, 6e-4, 4)), volumes, energies),
            BirchMurnaghan3rd(224.444565, 6.250619105057268e-4, 3.740369, -323.417714),
        )
        @test _isapprox(
            linfit(EnergyEOS(BirchMurnaghan3rd(224, 6e-4, 4)), volumes, energies),
            BirchMurnaghan3rd(224.444565, 6.250619105057268e-4, 3.740369, -323.417714),
        )
        # See https://github.com/aoterodelaroza/asturfit/blob/4de9b41/test/test03.out#L30-L36
        @test _isapprox(
            nonlinfit(
                EnergyEOS(BirchMurnaghan4th(224, 6e-4, 4, -5460)),  # bohr^3, Ry/bohr^3, 1, bohr^3/Ry, Ry
                volumes,
                energies,
            ),
            BirchMurnaghan4th(
                224.457562,  # bohr^3
                6.229381129795094e-4,  # Ry/bohr^3
                3.730992,
                -5322.7030547560435,  # bohr^3/Ry
                -323.417712,  # Ry
            );
            rtol = 1e-5,
        )
        # See https://github.com/aoterodelaroza/asturfit/blob/4de9b41/test/test03.out#L30-L36
        # @test _isapprox(
        #     linfit(
        #         EnergyEOS(BirchMurnaghan4th(224, 6e-4, 4, -5460)),  # bohr^3, Ry/bohr^3, 1, bohr^3/Ry, Ry
        #         volumes,
        #         energies,
        #     ),
        #     BirchMurnaghan4th(
        #         224.457562,  # bohr^3
        #         6.229381129795094e-4,  # Ry/bohr^3
        #         3.730992,
        #         -5322.7030547560435,  # bohr^3/Ry
        #         -323.417712,  # Ry
        #     );
        #     rtol = 1e-5,
        # )
        # See https://github.com/aoterodelaroza/asturfit/blob/4de9b41/test/test03.out#L98-L105
        # @test _isapprox(
        #     nonlinfit(
        # EnergyEOS(        #         BirchMurnaghan5th(224.445371, 0.0006, 4, -5500, 3.884535907971559e7, -323)),
        #         volumes,
        #         energies,
        #     ),
        #     BirchMurnaghan5th(
        #         224.451813,
        #         0.0006228893043314733,
        #         3.736723,
        #         -5292.414119096362,
        #         6.3542116050611705e7,
        #         -323.417712,
        #     );
        #     atol = 1,
        # )  # FIXME: result is wrong
        # # No reference data, I run on my computer.
        @test _isapprox(
            nonlinfit(EnergyEOS(Vinet(224, 6e-4, 4)), volumes, energies),
            Vinet(
                224.45278665796354,
                6.313500637481759e-4,
                3.7312381477678853,
                -323.4177229576912,
            ),
        )
        @test _isapprox(
            nonlinfit(EnergyEOS(PoirierTarantola3rd(224, 6e-4, 4)), volumes, energies),
            PoirierTarantola3rd(224.509208, 6.3589226415983795e-4, 3.690448, -323.41773);
            atol = 1e-5,
        )
        # # FIXME: This cannot go through
        # @test _isapprox(
        #     nonlinfit(
        # EnergyEOS(        #         PoirierTarantola4th(220, 0.0006, 3.7, -5500, -323)),
        #         volumes,
        #         energies,
        #     ),
        #     PoirierTarantola4th(
        #         224.430182,
        #         0.0006232241765069493,
        #         3.758360,
        #         -5493.859729817176,
        #         -323.417712,
        #     ),
        # )
        # See https://github.com/aoterodelaroza/asturfit/blob/4de9b41/test/test03.out#L98-L105
        # @test _isapprox(
        #     nonlinfit(
        #         EnergyEOS(PoirierTarantola5th(
        #             224.445371,
        #             0.0006,
        #             3.8,
        #             -5500,
        #             6e7,
        #             -323,
        #         )),
        #         volumes,
        #         energies,
        #     ),
        #     PoirierTarantola5th(
        #         224.451250,
        #         0.000622877204137392,
        #         3.737484,
        #         -5283.999708607125,
        #         6.296000262990379e7,
        #         -323.417712,
        #     );
        #     rtol = 0.05,
        # )
    end

    @testset "with units" begin
        volumes, energies = data["volume"] * u"bohr^3", data["energy"] * u"Ry"
        # See https://github.com/aoterodelaroza/asturfit/blob/4de9b41/test/test03.out#L117-L122
        @test _isapprox(
            nonlinfit(EnergyEOS(Murnaghan(224u"bohr^3", 9u"GPa", 4)), volumes, energies),
            Murnaghan(224.501825u"bohr^3", 8.896845u"GPa", 3.723835, -323.417686u"Ry");
            atol = 1e-5,
        )
        # See https://github.com/aoterodelaroza/asturfit/blob/4de9b41/test/test03.out#L15-L20
        @test _isapprox(
            nonlinfit(
                EnergyEOS(BirchMurnaghan3rd(224u"bohr^3", 9u"GPa", 4)),
                volumes,
                energies,
            ),
            BirchMurnaghan3rd(
                224.444565u"bohr^3",
                9.194978u"GPa",
                3.740369,
                -323.417714u"Ry",
            ),
        )
        # See https://github.com/aoterodelaroza/asturfit/blob/4de9b41/test/test03.out#L30-L36
        @test _isapprox(
            nonlinfit(
                EnergyEOS(BirchMurnaghan4th(224u"bohr^3", 9u"GPa", 4, -0.3u"1/GPa")),
                volumes,
                energies,
            ),
            BirchMurnaghan4th(
                224.457562u"bohr^3",
                9.163736u"GPa",
                3.730992,
                -0.361830u"1/GPa",
                -323.417712u"Ry",  # Ry
            );
            rtol = 1e-7,
        )
        # See https://github.com/aoterodelaroza/asturfit/blob/4de9b41/test/test03.out#L66-L71
        @test _isapprox(
            nonlinfit(
                EnergyEOS(PoirierTarantola3rd(224u"bohr^3", 9u"GPa", 4)),
                volumes,
                energies,
            ),
            PoirierTarantola3rd(
                224.509208u"bohr^3",
                9.354298u"GPa",
                3.690448,
                -323.41773u"Ry",
            );
            atol = 1e-5,
        )
    end
end

end
