module Collections

using IntervalArithmetic
using Measurements: Measurement, measurement
using Test: @test, @testset
using Unitful: Quantity, @u_str

using EquationsOfStateOfSolids.Collections

@testset "Test EOS promotion" begin
    @test typeof(Murnaghan(1, 2, 3.0, 0)) === Murnaghan{Float64}
    @test typeof(BirchMurnaghan2nd(1, 2.0, 0)) === BirchMurnaghan2nd{Float64}
    @test typeof(BirchMurnaghan3rd(1, 2, 3.0, 0)) === BirchMurnaghan3rd{Float64}
    @test typeof(BirchMurnaghan4th(1, 2.0, 3, 4, 0)) === BirchMurnaghan4th{Float64}
    @test typeof(PoirierTarantola2nd(1, 2.0, 0)) === PoirierTarantola2nd{Float64}
    @test typeof(PoirierTarantola3rd(1, 2, 3.0, 0)) === PoirierTarantola3rd{Float64}
    @test typeof(PoirierTarantola4th(1, 2, 3, 4, 0)) === PoirierTarantola4th{Int}
    @test typeof(Vinet(1, 2, 3.0, 0)) === Vinet{Float64}
    # @test typeof(AntonSchmidt(1, 2, 3.0, 0)) === AntonSchmidt{Float64}
    # @test typeof(BreenanStacey(1, 2, 3.0, 0)) === BreenanStacey{Float64}
    @test Murnaghan(1, Int32(2), Int8(3), 0) === Murnaghan{Int}(1, 2, 3, 0)
    @test Murnaghan(1, 2 // 1, Int8(3), 0) ===
          Murnaghan{Rational{Int}}(1 // 1, 2 // 1, 3 // 1, 0 // 1)
    @test BirchMurnaghan2nd(1, Int32(2), 0) === BirchMurnaghan2nd{Int}(1, 2, 0)
    @test BirchMurnaghan2nd(1 // 1, Int32(2), 0) ===
          BirchMurnaghan2nd{Rational{Int}}(1 // 1, 2 // 1, 0 // 1)
    @test BirchMurnaghan3rd(Int8(1), 2, 4, 0) === BirchMurnaghan3rd{Int}(1, 2, 4, 0)
    @test BirchMurnaghan3rd(Int8(1), 2 // 1, 4, 0) ===
          BirchMurnaghan3rd{Rational{Int}}(1 // 1, 2 // 1, 4 // 1, 0 // 1)
    @test BirchMurnaghan4th(Int8(1), 2, 4, Int16(5), 6) ===
          BirchMurnaghan4th{Int}(1, 2, 4, 5, 6)
    @test BirchMurnaghan4th(Int8(1), 2 // 1, 4, Int16(5), 6) ===
          BirchMurnaghan4th{Rational{Int}}(1 // 1, 2 // 1, 4 // 1, 5 // 1, 6 // 1)
    @test Murnaghan{Float32}(Int(1), 2 // 1, Int8(3), Float64(4)) ===
          Murnaghan(Float32(1.0), Float32(2.0), Float32(3.0), Float32(4.0))
    @test Murnaghan{BigFloat}(Int(1), 2 // 1, Int8(3), Float64(4)) ==
          Murnaghan(BigFloat(1.0), BigFloat(2.0), BigFloat(3.0), BigFloat(4.0))
    @test BirchMurnaghan4th(Int8(1), 2 // 1, big(4), Int16(5), 6) ==
          BirchMurnaghan4th{Rational{BigInt}}(1 // 1, 2 // 1, 4 // 1, 5 // 1, 6 // 1)
    @test BirchMurnaghan4th(Int8(1), 2 // 1, big(4.0), Int16(5), 6) ==
          BirchMurnaghan4th{BigFloat}(1.0, 2.0, 4.0, 5.0, 6.0)
    @test BirchMurnaghan4th(Int8(1), 2 // 1, big(4), Int16(5), 6.0) ==
          BirchMurnaghan4th{BigFloat}(1.0, 2.0, 4.0, 5.0, 6.0)
    @test PoirierTarantola2nd(Int8(1), 2, 3) === PoirierTarantola2nd{Int}(1, 2, 3)
    @test PoirierTarantola2nd(Int8(1), 2 // 1, 3) ===
          PoirierTarantola2nd{Rational{Int}}(1 // 1, 2 // 1, 3 // 1)
    @test PoirierTarantola3rd(Int8(1), 2, 3, Int16(4)) ===
          PoirierTarantola3rd{Int}(1, 2, 3, 4)
    @test PoirierTarantola3rd(Int8(1), 2 // 1, 3 // 1, Int16(4)) ===
          PoirierTarantola3rd{Rational{Int}}(1 // 1, 2 // 1, 3 // 1, 4 // 1)
    @test PoirierTarantola4th(Int8(1), 2, 3, Int16(4), 5) ===
          PoirierTarantola4th{Int}(1, 2, 3, 4, 5)
    @test PoirierTarantola4th(Int8(1), 2 // 1, 3, Int16(4), 5) ===
          PoirierTarantola4th{Rational{Int}}(1 // 1, 2 // 1, 3 // 1, 4 // 1, 5 // 1)
    @test Vinet(Int8(1), 2, 3, Int16(4)) === Vinet{Int}(1, 2, 3, 4)
    @test Vinet(Int8(1), 2 // 1, 3, Int16(4)) ===
          Vinet{Rational{Int}}(1 // 1, 2 // 1, 3 // 1, 4 // 1)
    # @test AntonSchmidt(Int8(1), 2, 3, 0) === AntonSchmidt{Int}(1, 2, 3, 0)
    # @test AntonSchmidt(Int8(1), 2 // 1, 3, 0) ===
    #       AntonSchmidt{Rational{Int}}(1 // 1, 2 // 1, 3 // 1, 0 // 1)
    # @test BreenanStacey(1, 2, 3, 0) === BreenanStacey{Int}(1, 2, 3, 0)
    # @test BreenanStacey(1, 2, 3 // 1, 0) ===
    #       BreenanStacey{Rational{Int}}(1 // 1, 2 // 1, 3 // 1, 0 // 1)
    @test typeof(Murnaghan(1 * u"angstrom^3", 2 * u"eV/angstrom^3", 3.0, 4 * u"eV")) ===
          Murnaghan{Quantity{Float64}}
    @test typeof(Murnaghan(1 * u"angstrom^3", 2 * u"eV/angstrom^3", 3 // 2, 4 * u"eV")) ===
          Murnaghan{Quantity{Rational{Int}}}
    @test typeof(BirchMurnaghan2nd(
        1 * u"angstrom^3",
        2 * u"eV/angstrom^3",
        3.0 * u"eV",
    )) === BirchMurnaghan2nd{Quantity{Float64}}
    @test typeof(BirchMurnaghan2nd(
        1 * u"angstrom^3",
        2 * u"eV/angstrom^3",
        3 // 1 * u"eV",
    )) === BirchMurnaghan2nd{Quantity{Rational{Int}}}
    @test typeof(BirchMurnaghan3rd(1 * u"angstrom^3", 2 * u"GPa", 4.0, 3 * u"eV")) ===
          BirchMurnaghan3rd{Quantity{Float64}}
    @test typeof(BirchMurnaghan3rd(1 * u"angstrom^3", 2 * u"GPa", 4 // 1, 3 * u"eV")) ===
          BirchMurnaghan3rd{Quantity{Rational{Int}}}
    @test typeof(BirchMurnaghan4th(
        1 * u"nm^3",
        2 * u"GPa",
        3.0,
        4 * u"1/GPa",
        5 * u"eV",
    )) === BirchMurnaghan4th{Quantity{Float64}}
    @test typeof(BirchMurnaghan4th(
        1 * u"nm^3",
        2 * u"GPa",
        3 // 1,
        4 * u"1/GPa",
        5 * u"eV",
    )) === BirchMurnaghan4th{Quantity{Rational{Int}}}
    @test typeof(PoirierTarantola2nd(1 * u"nm^3", 2 * u"GPa", 3.0 * u"eV")) ===
          PoirierTarantola2nd{Quantity{Float64}}
    @test typeof(PoirierTarantola2nd(1 * u"nm^3", 2 * u"GPa", 3 // 1 * u"eV")) ===
          PoirierTarantola2nd{Quantity{Rational{Int}}}
    @test typeof(PoirierTarantola3rd(1 * u"nm^3", 2 * u"GPa", 3, 4.0 * u"eV")) ===
          PoirierTarantola3rd{Quantity{Float64}}
    @test typeof(PoirierTarantola3rd(1 * u"nm^3", 2 * u"GPa", 3, 4 // 1 * u"eV")) ===
          PoirierTarantola3rd{Quantity{Rational{Int}}}
    @test typeof(PoirierTarantola4th(
        1 * u"nm^3",
        2 * u"GPa",
        3,
        4 * u"1/GPa",
        5.0 * u"eV",
    )) === PoirierTarantola4th{Quantity{Float64}}
    @test typeof(PoirierTarantola4th(
        1 * u"nm^3",
        2 * u"GPa",
        3,
        4 * u"1/GPa",
        5 // 1 * u"eV",
    )) === PoirierTarantola4th{Quantity{Rational{Int}}}
    @test typeof(Vinet(1 * u"nm^3", 2 * u"GPa", 3, 4.0 * u"eV")) ===
          Vinet{Quantity{Float64}}
    @test typeof(Vinet(1 * u"nm^3", 2 * u"GPa", 3, 4 // 1 * u"eV")) ===
          Vinet{Quantity{Rational{Int}}}
    # @test typeof(AntonSchmidt(1 * u"nm^3", 2 * u"GPa", 3.0, 0 * u"eV")) ===
    #       AntonSchmidt{Quantity{Float64}}
    # @test typeof(AntonSchmidt(1 * u"nm^3", 2 * u"GPa", 3 // 1, 0 * u"eV")) ===
    #       AntonSchmidt{Quantity{Rational{Int}}}
    # @test typeof(BreenanStacey(1 * u"nm^3", 2 * u"GPa", 3.0, 0 * u"eV")) ===
    #       BreenanStacey{Quantity{Float64}}
    # @test typeof(BreenanStacey(1 * u"nm^3", 2 * u"GPa", 3 // 1, 0 * u"eV")) ===
    #   BreenanStacey{Quantity{Rational{Int}}}
    @test BirchMurnaghan3rd(1 * u"angstrom^3", 2 * u"GPa", 4 // 1, 3 * u"eV").b′0 isa
          DimensionlessQuantity
end

@testset "Test default EOS parameter `e0` and promotion" begin
    @test Murnaghan(1, 2, 3.0) === Murnaghan(1.0, 2.0, 3.0, 0.0)
    @test BirchMurnaghan2nd(1, 2.0) === BirchMurnaghan2nd(1.0, 2.0, 0.0)
    @test BirchMurnaghan3rd(1, 2, 3.0) === BirchMurnaghan3rd(1.0, 2.0, 3.0, 0.0)
    @test BirchMurnaghan4th(1, 2.0, 3, 4) === BirchMurnaghan4th(1.0, 2.0, 3.0, 4.0, 0.0)
    @test Vinet(1, 2, 3.0) === Vinet(1.0, 2.0, 3.0, 0.0)
    @test PoirierTarantola2nd(1, 2.0) === PoirierTarantola2nd(1.0, 2.0, 0.0)
    @test PoirierTarantola3rd(1, 2, 3.0) === PoirierTarantola3rd(1.0, 2.0, 3.0, 0.0)
    @test PoirierTarantola4th(1, 2, 3, 4) === PoirierTarantola4th(1, 2, 3, 4, 0)
    # @test AntonSchmidt(1, 2, 3.0) === AntonSchmidt(1.0, 2.0, 3.0, 0.0)
    # @test BreenanStacey(1, 2, 3.0) === BreenanStacey(1.0, 2.0, 3.0, 0.0)
    @test typeof(Murnaghan(1 * u"angstrom^3", 2 * u"eV/angstrom^3", 3)) ===
          Murnaghan{Quantity{Int}}
    @test typeof(Murnaghan(1 * u"angstrom^3", 2 * u"eV/angstrom^3", 3.0)) ===
          Murnaghan{Quantity{Float64}}
    @test typeof(BirchMurnaghan2nd(1 * u"nm^3", 2 * u"GPa")) ===
          BirchMurnaghan2nd{Quantity{Int}}
    @test typeof(BirchMurnaghan2nd(1 * u"nm^3", 2.0 * u"GPa")) ===
          BirchMurnaghan2nd{Quantity{Float64}}
    @test typeof(BirchMurnaghan3rd(1 * u"nm^3", 2 * u"GPa", 4)) ===
          BirchMurnaghan3rd{Quantity{Int}}
    @test typeof(BirchMurnaghan3rd(1 * u"nm^3", 2 * u"GPa", 4.0)) ===
          BirchMurnaghan3rd{Quantity{Float64}}
    @test typeof(BirchMurnaghan4th(1 * u"nm^3", 2 * u"GPa", 3, 4 * u"1/GPa")) ===
          BirchMurnaghan4th{Quantity{Int}}
    @test typeof(BirchMurnaghan4th(1 * u"nm^3", 2 * u"GPa", 3, 4.0 * u"1/GPa")) ===
          BirchMurnaghan4th{Quantity{Float64}}
    @test typeof(PoirierTarantola2nd(1 * u"nm^3", 2 * u"GPa")) ===
          PoirierTarantola2nd{Quantity{Int}}
    @test typeof(PoirierTarantola2nd(1 * u"nm^3", 2.0 * u"GPa")) ===
          PoirierTarantola2nd{Quantity{Float64}}
    @test typeof(PoirierTarantola3rd(1 * u"nm^3", 2 * u"GPa", 3)) ===
          PoirierTarantola3rd{Quantity{Int}}
    @test typeof(PoirierTarantola3rd(1 * u"nm^3", 2 * u"GPa", 3.0)) ===
          PoirierTarantola3rd{Quantity{Float64}}
    @test typeof(PoirierTarantola4th(1 * u"nm^3", 2 * u"GPa", 3, 4 * u"1/GPa")) ===
          PoirierTarantola4th{Quantity{Int}}
    @test typeof(PoirierTarantola4th(1 * u"nm^3", 2 * u"GPa", 3, 4.0 * u"1/GPa")) ===
          PoirierTarantola4th{Quantity{Float64}}
    @test typeof(Vinet(1 * u"nm^3", 2 * u"GPa", 3)) === Vinet{Quantity{Int}}
    @test typeof(Vinet(1 * u"nm^3", 2 * u"GPa", 3.0)) === Vinet{Quantity{Float64}}
    #@test typeof(AntonSchmidt(1u"nm^3", 2u"GPa", 3))
    #@test typeof(AntonSchmidt(1u"nm^3", 2u"GPa", 3.0))
    #@test typeof(BreenanStacey(1u"nm^3", 2u"GPa", 3))
    #@test typeof(BreenanStacey(1u"nm^3", 2u"GPa", 3.0)
    # @test typeof(PolynomialEOS(1, [1, 2, 3], 4.0)) === PolynomialEOS{3,Float64}
    # @test typeof(PolynomialEOS(1, [1.0, 2, 3, 4], 4)) === PolynomialEOS{4,Float64}
    # @test typeof(PolynomialEOS(1.0, [1, 2, 3, 4])) === PolynomialEOS{4,Float64}
    # @test typeof(PolynomialEOS(
    #     1 * u"nm^3",
    #     [1 * u"eV/nm^3", 2, 3 * u"nm^3/eV"],
    #     4 * u"eV",
    # )) === PolynomialEOS{3,Quantity{Int}}
    # @test typeof(PolynomialEOS(
    #     1 * u"nm^3",
    #     [1 * u"eV/nm^3", 2, 3 * u"nm^3/eV"],
    #     4.0 * u"eV",
    # )) === PolynomialEOS{3,Quantity{Float64}}
end

@testset "`float` on an EOS" begin
    @test float(Vinet(1, 2, 3)) == Vinet(1.0, 2.0, 3.0, 0.0)
    @test float(Vinet(1 * u"nm^3", 2 * u"GPa", 3)) ==
          Vinet(1.0 * u"nm^3", 2.0 * u"GPa", 3.0)
    # @test float(PolynomialEOS(1, [1, 2, 3])) == PolynomialEOS(1.0, [1.0, 2.0, 3.0])
    # @test float(PolynomialEOS(1 * u"nm^3", [1 * u"eV/nm^3", 2, 3 * u"nm^3/eV"])) ==
    #       PolynomialEOS(1.0 * u"nm^3", [1.0 * u"eV/nm^3", 2.0, 3.0 * u"nm^3/eV"])
end # testset

@testset "Test counstruction" begin
    @test typeof(BirchMurnaghan(1u"angstrom^3", 3u"GPa", 4.0)) ==
          BirchMurnaghan{3,Quantity{Float64}}
    @test typeof(BirchMurnaghan(
        measurement("1 +- 0.1"),
        3 // 1,
        2,
        measurement("-123.4(56)"),
    )) == BirchMurnaghan{4,Measurement{Float64}}
    @test typeof(BirchMurnaghan(1..3, 2, 2..4)) == BirchMurnaghan{3,Interval{Float64}}
end

# Data in the following tests are from
# https://github.com/materialsproject/pymatgen/blob/1f0957b8525ddc7d12ea348a19caecebe6c7ff34/pymatgen/analysis/tests/test_eos.py
@testset "Test data from Pymatgen" begin
    volumes = [
        25.987454833,
        26.9045702104,
        27.8430241908,
        28.8029649591,
        29.7848370694,
        30.7887887064,
        31.814968055,
        32.8638196693,
        33.9353435494,
        35.0299842495,
        36.1477417695,
        37.2892088485,
        38.4543854865,
        39.6437162376,
        40.857201102,
        42.095136449,
        43.3579668329,
        44.6456922537,
        45.9587572656,
        47.2973100535,
        48.6614988019,
        50.0517680652,
        51.4682660281,
        52.9112890601,
        54.3808371612,
        55.8775030703,
        57.4014349722,
        58.9526328669,
    ]
    energies = [
        -7.63622156576,
        -8.16831294894,
        -8.63871612686,
        -9.05181213218,
        -9.41170988374,
        -9.72238224345,
        -9.98744832526,
        -10.210309552,
        -10.3943401353,
        -10.5427238068,
        -10.6584266073,
        -10.7442240979,
        -10.8027285713,
        -10.8363890521,
        -10.8474912964,
        -10.838157792,
        -10.8103477586,
        -10.7659387815,
        -10.7066179666,
        -10.6339907853,
        -10.5495538639,
        -10.4546677714,
        -10.3506386542,
        -10.2386366017,
        -10.1197772808,
        -9.99504030111,
        -9.86535084973,
        -9.73155247952,
    ]
    x = BirchMurnaghan3rd(40.98926572528106, 0.5369258245417454, 4.178644235500821)
end

end
