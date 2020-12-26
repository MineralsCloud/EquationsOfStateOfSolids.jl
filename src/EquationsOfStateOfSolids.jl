module EquationsOfStateOfSolids

_ispositive(x) = x > zero(x)  # Do not export!

include("collections.jl")
include("Fitting.jl")
include("Volume.jl")

end
