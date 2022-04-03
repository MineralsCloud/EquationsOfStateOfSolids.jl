module FiniteStrains

using Test: @test, @testset, @test_throws

using EquationsOfStateOfSolids.FiniteStrains: ToEulerianStrain

# See issue https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/issues/146
@testset "Test `ToEulerianStrain`" begin
    volumes =
        ([17.789658, 18.382125, 18.987603, 19.336585, 19.606232, 20.238155, 20.883512])
    v0 = 19.361802843683524
    strains = map(ToEulerianStrain(v0), volumes)
    @test first(strains) == 0.029040354449144212
end

end
