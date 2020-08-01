module Fitting

using ConstructionBase: constructorof
using LsqFit: curve_fit
using Setfield
using Unitful: AbstractQuantity, NoDims, upreferred, ustrip, unit, dimension, @u_str

using ..Collections

export linfit, nonlinfit

function nonlinfit(eos, xs, ys; kwargs...)
    params = [eos.v0, eos.b0, eos.b′0]
    function model(x, p)
        @set! eos.v0 = p[1]
        @set! eos.b0 = p[2]
        @set! eos.b′0 = p[3]
        return eos.(x)
    end
    fit = curve_fit(model, xs, ys, params; kwargs...)
    return fit
end # function nonlinfit

end
