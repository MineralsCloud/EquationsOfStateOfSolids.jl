using EquationsOfStateOfSolids
using Test

@testset "EquationsOfStateOfSolids.jl" begin
    include("collections.jl")
    include("FiniteStrains.jl")
    include("Fitting.jl")
end
