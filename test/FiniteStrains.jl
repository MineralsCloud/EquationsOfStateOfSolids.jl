module FiniteStrains

using Test: @test, @testset, @test_throws

using EquationsOfStateOfSolids.FiniteStrains: ToEulerianStrain

# See issue https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/issues/146
@testset "Test `ToEulerianStrain`" begin
    v = 17.789658
    v0 = 19.361802843683524
    @test ToEulerianStrain(v0)(v) == 0.029040354449144323
end

end
