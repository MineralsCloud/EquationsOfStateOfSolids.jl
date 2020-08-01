module Fitting

using ConstructionBase: constructorof
using LsqFit: curve_fit
using LeastSquaresOptim
using Setfield
using Unitful: AbstractQuantity, NoDims, upreferred, ustrip, unit, dimension, @u_str

using ..Collections

export linfit, nonlinfit

function nonlinfit(eos, xs, ys; kwargs...)
    params = [eos.v0, eos.b0, eos.b′0, 0]
    f = costfunction(xs, ys)
    return optimize(f, params, LevenbergMarquardt())
end # function nonlinfit

end
