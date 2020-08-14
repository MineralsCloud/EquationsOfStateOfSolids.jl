module Fitting

using ConstructionBase: constructorof
using LsqFit: curve_fit
using Unitful: AbstractQuantity, NoDims, upreferred, ustrip, unit, dimension, @u_str

using ..Collections: EnergyEoss, PressureEoss, fieldvalues

export linfit, nonlinfit

const float_collect = float ∘ collect

function nonlinfit(eos::S, xs, ys; kwargs...) where {T,S<:EnergyEoss{T}}
    constructor = constructorof(S) ∘ constructorof(T)
    # model(x, p) = map(constructor(p), x)
    model(x, p) = [constructor(p[1:end-1])(xx, p[end]) for xx in x]
    fit = curve_fit(
        model,
        float_collect(xs),
        float_collect(ys),
        vcat(float_collect(fieldvalues(eos.params)), minimum(ys));
        kwargs...,
    )
    if fit.converged
        return constructorof(T)(fit.param), fit
    else
        @error "fitting is not converged! Check your initial parameters!"
        return nothing, fit
    end
end # function nonlinfit
function nonlinfit(eos::S, xs, ys; kwargs...) where {T,S<:PressureEoss{T}}
    constructor = constructorof(S) ∘ constructorof(T)
    model(x, p) = map(constructor(p), x)
    fit = curve_fit(
        model,
        float_collect(xs),
        float_collect(ys),
        float_collect(fieldvalues(eos.params));
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
