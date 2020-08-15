module Fitting

using ConstructionBase: constructorof
using LsqFit: curve_fit
using Unitful: AbstractQuantity, NoDims, upreferred, ustrip, unit, dimension, @u_str

using ..Collections: EquationOfStateOfSolids

export linfit, nonlinfit

function nonlinfit(eos::EquationOfStateOfSolids{T}, xs, ys; kwargs...) where {T}
    model = _create_model(eos)
    fit =
        curve_fit(model, float.(xs), float.(ys), float.(fieldvalues(eos.param)); kwargs...)
    if fit.converged
        return constructorof(T)(fit.param), fit
    else
        @error "fitting is not converged! Check your initial parameters!"
        return nothing, fit
    end
end # function nonlinfit

function _create_model(::S) where {T,S<:EquationOfStateOfSolids{T}}
    constructor = constructorof(S) ∘ constructorof(T)
    return (x, p) -> map(constructor(p), x)
end

fieldvalues(x) = (getfield(x, i) for i in 1:nfields(x))

end
