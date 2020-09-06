module Collections

using Measurements: Measurement, measurement
using SymPy: Sym, symbols
using Test: @test, @testset
using Unitful: Quantity, DimensionlessQuantity, @u_str

using EquationsOfStateOfSolids.Collections

@testset "Promoting `eltype`" begin
    @testset "Promoting to floating-point numbers" begin
        @test eltype(Murnaghan(1, 2, 3.0, 0)) === Float64
        @test eltype(BirchMurnaghan2nd(1, 2.0, 0)) === Float64
        @test eltype(BirchMurnaghan3rd(1, 2, 3.0, 0)) === Float64
        @test eltype(BirchMurnaghan4th(1, 2.0, 3, 4, 0)) === Float64
        @test eltype(PoirierTarantola2nd(1, 2.0, 0)) === Float64
        @test eltype(PoirierTarantola3rd(1, 2, 3.0, 0)) === Float64
        @test eltype(Vinet(1, 2, 3.0, 0)) === Float64
        @test eltype(AntonSchmidt(1, 2, 3.0, 0)) === Float64
        # @test eltype(BreenanStacey(1, 2, 3.0, 0)) === Float64
        @test eltype(Murnaghan{Float32}(Int(1), 2 // 1, Int8(3), Float64(4))) === Float32
        @test eltype(Murnaghan{BigFloat}(Int(1), 2 // 1, Int8(3), Float64(4))) === BigFloat
        @test eltype(PoirierTarantola4th{Float16}(Int8(1), 2 // 1, 4, Int16(5), 6)) ===
              Float16
        @test eltype(BirchMurnaghan4th(Int8(1), 2 // 1, big(4.0), Int16(5), 6)) === BigFloat
        @test eltype(BirchMurnaghan4th(Int8(1), 2 // 1, big(4), Int16(5), 6.0)) === BigFloat
    end

    @testset "Promoting to rationals" begin
        @test eltype(PoirierTarantola4th(1, 2, 3, 4, 0)) === PoirierTarantola4th{Int}
        @test Murnaghan(Int32(1), Int16(2), Int8(3), 0) === Murnaghan{Int}(1, 2, 3, 0)
        @test Murnaghan(1, 2 // 1, Int8(3), 0) ===
              Murnaghan{Rational{Int}}(1 // 1, 2 // 1, 3 // 1, 0 // 1)
        @test BirchMurnaghan2nd(1, Int8(2), 0) === BirchMurnaghan2nd{Int}(1, 2, 0)
        @test BirchMurnaghan2nd(1 // 1, Int32(2)) ===
              BirchMurnaghan2nd{Rational{Int}}(1 // 1, 2 // 1, 0 // 1)
        @test BirchMurnaghan3rd(Int8(1), 2, 4, 0) === BirchMurnaghan3rd{Int}(1, 2, 4, 0)
        @test BirchMurnaghan3rd(Int8(1), 2 // 1, 4, 0) ===
              BirchMurnaghan3rd{Rational{Int}}(1 // 1, 2 // 1, 4 // 1, 0 // 1)
        @test BirchMurnaghan4th(Int8(1), 2, 4, Int16(5), 6) ===
              BirchMurnaghan4th{Int}(1, 2, 4, 5, 6)
        @test BirchMurnaghan4th(Int8(1), 2 // 1, 4, Int16(5), 6) ===
              BirchMurnaghan4th{Rational{Int}}(1 // 1, 2 // 1, 4 // 1, 5 // 1, 6 // 1)
        @test BirchMurnaghan4th(Int8(1), 2 // 1, big(4), Int16(5), 6) ==
              BirchMurnaghan4th{Rational{BigInt}}(1 // 1, 2 // 1, 4 // 1, 5 // 1, 6 // 1)
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
        @test AntonSchmidt(Int8(1), 2, 3, 0) === AntonSchmidt{Int}(1, 2, 3, 0)
        @test AntonSchmidt(Int8(1), 2 // 1, 3, 0) ===
              AntonSchmidt{Rational{Int}}(1 // 1, 2 // 1, 3 // 1, 0 // 1)
        #   @test BreenanStacey(1, 2, 3, 0) === BreenanStacey{Int}(1, 2, 3, 0)
        #   @test BreenanStacey(1, 2, 3 // 1, 0) ===
        #         BreenanStacey{Rational{Int}}(1 // 1, 2 // 1, 3 // 1, 0 // 1)
    end

    @testset "Promoting with units" begin
        @test Murnaghan(1u"angstrom^3", 2u"eV/angstrom^3", 3.0, 4u"eV") ===
              Murnaghan(1.0u"angstrom^3", 2.0u"eV/angstrom^3", 3.0, 4.0u"eV")
        @test Murnaghan(1u"angstrom^3", 2u"eV/nm^3", 3 // 2, 4u"eV") ===
              Murnaghan((1 // 1)u"angstrom^3", (2 // 1)u"eV/nm^3", 3 // 2, (4 // 1)u"eV")
        @test BirchMurnaghan2nd(1u"angstrom^3", 2u"eV/angstrom^3", 3.0u"J") ===
              BirchMurnaghan2nd(1.0u"angstrom^3", 2.0u"eV/angstrom^3", 3.0u"J")
        @test BirchMurnaghan2nd(1u"pm^3", 2u"eV/angstrom^3", (3 // 1)u"eV") ===
              BirchMurnaghan2nd((1 // 1)u"pm^3", (2 // 1)u"eV/angstrom^3", (3 // 1)u"eV")
        @test BirchMurnaghan3rd(1u"angstrom^3", 2u"GPa", 4.0, 3u"eV") ===
              BirchMurnaghan3rd(1.0u"angstrom^3", 2.0u"GPa", 4.0, 3.0u"eV")
        @test BirchMurnaghan3rd(1u"angstrom^3", 2u"GPa", 4 // 1, 3u"eV") ===
              BirchMurnaghan3rd(
            (1 // 1)u"angstrom^3",
            (2 // 1)u"GPa",
            4 // 1,
            (3 // 1)u"eV",
        )
        @test BirchMurnaghan4th(1u"nm^3", 2u"GPa", 3.0, 4u"GPa^-1", 5u"eV") ===
              BirchMurnaghan4th(1.0u"nm^3", 2.0u"GPa", 3.0, 4.0u"GPa^-1", 5.0u"eV")
        @test BirchMurnaghan4th(1u"nm^3", 2u"GPa", 3 // 1, 4u"1/GPa", 5u"J") ===
              BirchMurnaghan4th(
            (1 // 1)u"nm^3",
            (2 // 1)u"GPa",
            3 // 1,
            (4 // 1)u"1/GPa",
            (5 // 1)u"J",
        )
        @test PoirierTarantola2nd(1u"pm^3", 2u"GPa", 3.0u"eV") ===
              PoirierTarantola2nd(1.0u"pm^3", 2.0u"GPa", 3.0u"eV")
        @test PoirierTarantola2nd(1u"nm^3", 2u"GPa", (3 // 1)u"eV") ===
              PoirierTarantola2nd((1 // 1)u"nm^3", (2 // 1)u"GPa", (3 // 1)u"eV")
        @test PoirierTarantola3rd(1u"nm^3", 2u"GPa", 3, 4.0u"eV") ===
              PoirierTarantola3rd(1.0u"nm^3", 2.0u"GPa", 3, 4.0u"eV")
        @test PoirierTarantola3rd(1u"nm^3", 2u"GPa", 3, (4 // 1)u"eV") ===
              PoirierTarantola3rd((1 // 1)u"nm^3", (2 // 1)u"GPa", 3 // 1, (4 // 1)u"eV")
        @test PoirierTarantola4th(1u"nm^3", 2u"GPa", 3, 1 / 0.25u"Pa", 5.0u"eV") ===
              PoirierTarantola4th(1.0u"nm^3", 2.0u"GPa", 3.0, 4.0u"1/Pa", 5.0u"eV")
        @test PoirierTarantola4th(1u"nm^3", 2u"GPa", 3 // 1, 4u"GPa^(-1)", 5u"eV") ===
              PoirierTarantola4th(
            (1 // 1)u"nm^3",
            (2 // 1)u"GPa",
            3 // 1,
            (4 // 1)u"GPa^(-1)",
            (5 // 1)u"eV",
        )
        @test Vinet(1u"nm^3", 2u"GPa", 3, 4.0u"eV") ===
              Vinet(1.0u"nm^3", 2.0u"GPa", 3.0, 4.0u"eV")
        @test Vinet(1u"nm^3", 2u"GPa", 3, (4 // 1)u"eV") ===
              Vinet((1 // 1)u"nm^3", (2 // 1)u"GPa", 3 // 1, (4 // 1)u"eV")
        @test AntonSchmidt(1u"nm^3", 2u"GPa", 3.0, 4u"eV") ===
              AntonSchmidt(1.0u"nm^3", 2.0u"GPa", 3.0, 4.0u"eV")
        @test AntonSchmidt(1u"nm^3", 2u"GPa", 3 // 1, 4u"eV") ===
              AntonSchmidt((1 // 1)u"nm^3", (2 // 1)u"GPa", 3 // 1, (4 // 1)u"eV")
        # @test BreenanStacey(1u"nm^3", 2u"GPa", 3.0, 0u"eV") ===
        #       BreenanStacey{Quantity{Float64}}
        # @test BreenanStacey(1u"nm^3", 2u"GPa", 3 // 1, 0u"eV") ===
        #   BreenanStacey{Quantity{Rational{Int}}}
        @test BirchMurnaghan3rd(1u"angstrom^3", 2u"GPa", 4 // 1, 3u"eV").b′0 isa
              DimensionlessQuantity
        @test AntonSchmidt(1u"nm^3", 2u"GPa", 3.0, 4u"eV").b′0 isa DimensionlessQuantity
    end
end

@testset "Parameter `e0` and promotion" begin
    @test Murnaghan(1, 2, 3.0) === Murnaghan(1.0, 2.0, 3.0, 0.0)
    @test BirchMurnaghan2nd(1, 2.0) === BirchMurnaghan2nd(1.0, 2.0, 0.0)
    @test BirchMurnaghan3rd(1, 2, 3.0) === BirchMurnaghan3rd(1.0, 2.0, 3.0, 0.0)
    @test BirchMurnaghan4th(1, 2.0, 3, 4) === BirchMurnaghan4th(1.0, 2.0, 3.0, 4.0, 0.0)
    @test Vinet(1, 2, 3.0) === Vinet(1.0, 2.0, 3.0, 0.0)
    @test PoirierTarantola2nd(1, 2.0) === PoirierTarantola2nd(1.0, 2.0, 0.0)
    @test PoirierTarantola3rd(1, 2, 3.0) === PoirierTarantola3rd(1.0, 2.0, 3.0, 0.0)
    @test PoirierTarantola4th(1, 2, 3, 4) === PoirierTarantola4th(1, 2, 3, 4, 0)
    @test AntonSchmidt(1, 2, 3.0) === AntonSchmidt(1.0, 2.0, 3.0, 0.0)
    # @test BreenanStacey(1, 2, 3.0) === BreenanStacey(1.0, 2.0, 3.0, 0.0)
    @test eltype(Murnaghan(1u"angstrom^3", 2u"eV/angstrom^3", 3)) === Quantity{Int}
    @test eltype(Murnaghan(1u"angstrom^3", 2u"eV/angstrom^3", 3.0)) === Quantity{Float64}
    @test eltype(BirchMurnaghan2nd(1u"nm^3", 2u"GPa")) === Quantity{Int}
    @test eltype(BirchMurnaghan2nd(1u"nm^3", 2.0u"GPa")) === Quantity{Float64}
    @test eltype(BirchMurnaghan3rd(1u"nm^3", 2u"GPa", 4)) === Quantity{Int}
    @test eltype(BirchMurnaghan3rd(1u"nm^3", 2u"GPa", 4.0)) === Quantity{Float64}
    @test eltype(BirchMurnaghan4th(1u"nm^3", 2u"GPa", 3, 4u"1/GPa")) === Quantity{Int}
    @test eltype(BirchMurnaghan4th(1u"nm^3", 2u"GPa", 3, 4.0u"1/GPa")) === Quantity{Float64}
    @test eltype(PoirierTarantola2nd(1u"nm^3", 2u"GPa")) === Quantity{Int}
    @test eltype(PoirierTarantola2nd(1u"nm^3", 2.0u"GPa")) === Quantity{Float64}
    @test eltype(PoirierTarantola3rd(1u"nm^3", 2u"GPa", 3)) === Quantity{Int}
    @test eltype(PoirierTarantola3rd(1u"nm^3", 2u"GPa", 3.0)) === Quantity{Float64}
    @test eltype(PoirierTarantola4th(1u"nm^3", 2u"GPa", 3, 4u"1/GPa")) === Quantity{Int}
    @test eltype(PoirierTarantola4th(1u"nm^3", 2u"GPa", 3, 4.0u"1/GPa")) ===
          Quantity{Float64}
    @test eltype(Vinet(1u"nm^3", 2u"GPa", 3)) === Quantity{Int}
    @test eltype(Vinet(1u"nm^3", 2u"GPa", 3.0)) === Quantity{Float64}
    @test eltype(AntonSchmidt(1u"nm^3", 2u"GPa", 3)) === Quantity{Int}
    @test eltype(AntonSchmidt(1u"nm^3", 2u"GPa", 3.0)) === Quantity{Float64}
    #@test eltype(BreenanStacey(1u"nm^3", 2u"GPa", 3))
    #@test eltype(BreenanStacey(1u"nm^3", 2u"GPa", 3.0)
end

@testset "`float` on an EOS" begin
    @test float(Vinet(1, 2, 3)) == Vinet(1.0, 2.0, 3.0, 0.0)
    @test float(Vinet(1u"nm^3", 2u"GPa", 3)) == Vinet(1.0u"nm^3", 2.0u"GPa", 3.0)
end

@testset "Other element types" begin
    @test eltype(BirchMurnaghan4th(
        measurement("1 +- 0.1"),
        3 // 1,
        2,
        measurement("-123.4(56)"),
    )) === Measurement{Float64}
    v0, b0, b′0, b″0, e0 = symbols("v0, b0, b′0, b″0, e0")
    @test eltype(BirchMurnaghan4th(v0, b0, b′0, b″0, e0)) === Sym
end

end
