module Fitting

using ConstructionBase: constructorof
using LsqFit: curve_fit
using Unitful: AbstractQuantity, NoDims, upreferred, ustrip, unit, dimension, @u_str

using ..Collections

export linfit, nonlinfit

function nonlinfit(eos, xs, ys; kwargs...)
end # function nonlinfit

end
