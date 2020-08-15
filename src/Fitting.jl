module Fitting

using ConstructionBase: constructorof
using LsqFit: curve_fit
using Unitful: AbstractQuantity, NoDims, upreferred, ustrip, unit, dimension, @u_str

using ..Collections: EnergyEoss, PressureEoss

export linfit, nonlinfit

function nonlinfit(eos::S, xs, ys; kwargs...) where {T,S<:EnergyEoss{T}}
    constructor = constructorof(S) ∘ constructorof(T)
    model(x, p) = map(constructor(p), x)
    fit =
        curve_fit(model, float.(xs), float.(ys), float.(fieldvalues(eos.param)); kwargs...)
    if fit.converged
        return constructorof(T)(fit.param), fit
    else
        @error "fitting is not converged! Check your initial parameters!"
        return nothing, fit
    end
end # function nonlinfit


fieldvalues(x) = (getfield(x, i) for i in 1:nfields(x))

end
