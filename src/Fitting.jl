module Fitting

using ConstructionBase: constructorof
using LsqFit: curve_fit, coef
using Unitful: AbstractQuantity, NoDims, upreferred, ustrip, unit, dimension, @u_str

using ..Collections: EquationOfStateOfSolids, FiniteStrainEossParam, PressureEoss

export linfit, nonlinfit

function nonlinfit(eos::EquationOfStateOfSolids{T}, xs, ys; kwargs...) where {T}
    model = createmodel(eos)
    fit =
        curve_fit(model, float.(xs), float.(ys), float.(fieldvalues(eos.param)); kwargs...)
    if fit.converged
        param = constructorof(T)(coef(fit))
        checkparam(param)
        return param, fit
    else
        @error "fitting is not converged! Check your initial parameters!"
        return nothing, fit
    end
end # function nonlinfit

function createmodel(::S) where {T,S<:EquationOfStateOfSolids{T}}  # Do not export!
    constructor = constructorof(S) ∘ constructorof(T)
    return (x, p) -> map(constructor(p), x)
end

function checkparam(param::FiniteStrainEossParam)  # Do not export!
    if param.v0 <= zero(param.v0) || param.b0 <= zero(param.b0)
        @error "fitted `v0` or `b0` is negative!"
    end
    # if PressureEoss(param)(minimum(v)) >= param.b0
    #     @warn "use higher order EOS!"
    # end
end

fieldvalues(x) = (getfield(x, i) for i in 1:nfields(x))  # Do not export!

end
