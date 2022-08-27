module Fitting

using ConstructionBase: setproperties
using Unitful: ustrip, unit

using ..EquationsOfStateOfSolids: EnergyEquation

export fiteos

# See https://github.com/JuliaMath/Roots.jl/blob/bf0da62/src/utils.jl#L9-L11
struct ConvergenceFailed
    msg::String
end

abstract type FittingMethod end

include("linear.jl")
include("nonlinear.jl")

end
