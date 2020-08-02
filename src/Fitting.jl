module Fitting

using ConstructionBase: constructorof
using LsqFit: curve_fit
using Unitful: AbstractQuantity, NoDims, upreferred, ustrip, unit, dimension, @u_str

using ..Collections: EquationOfStateOfSolids

export linfit, nonlinfit

const float_collect = float ∘ collect

function nonlinfit(eos::S, xs, ys; kwargs...) where {T,S<:EquationOfStateOfSolids{T}}
    constructor = constructorof(S) ∘ constructorof(T)
    model(x, p) = map(constructor(p), x)
    fit = curve_fit(
        model,
        float_collect(xs),
        float_collect(ys),
        float_collect(eos.params.x0);
        kwargs...,
    )
    if fit.converged
        return constructorof(T)(fit.param), fit
    else
        @error "fitting is not converged! Check your initial parameters!"
        return nothing, fit
    end
end # function nonlinfit

end
