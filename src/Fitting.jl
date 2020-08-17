module Fitting

using ConstructionBase: constructorof
using LsqFit: curve_fit, coef
using Serialization: serialize
using Setfield: @set!

using ..Collections:
    EquationOfStateOfSolids, FiniteStrainEossParam, PressureEoss, EnergyEoss

export linfit, nonlinfit

function nonlinfit(
    eos::EquationOfStateOfSolids{T},
    xs,
    ys;
    xtol,
    gtol,
    maxiter::Integer = 1000,
    min_step_quality = 1e-3,
    good_step_quality = 0.75,
    silent = true,
    saveto = "",
) where {T}
    model = createmodel(eos)
    p0 = initparam(eos, ys)
    fit = curve_fit(
        model,
        float.(xs),
        float.(ys),
        float.(p0);
        x_tol = xtol,
        g_tol = gtol,
        maxIter = maxiter,
        min_step_quality = min_step_quality,
        good_step_quality = good_step_quality,
        show_trace = !silent,
    )
    if !isempty(saveto)
        open(saveto, "r+") do io
            serialize(io, fit)
        end
    end
    if fit.converged
        param = constructorof(T)(coef(fit))
        checkparam(param)
        return param
    else
        @error "fitting is not converged! Change initial parameters!"
        return nothing
    end
end # function nonlinfit

function createmodel(::S) where {T,S<:EquationOfStateOfSolids{T}}  # Do not export!
    constructor = constructorof(S) ∘ constructorof(T)
    return (x, p) -> map(constructor(p), x)
end

function checkparam(param::FiniteStrainParameters)  # Do not export!
    if param.v0 <= zero(param.v0) || param.b0 <= zero(param.b0)
        @error "fitted `v0` or `b0` is negative!"
    end
    # if PressureEoss(param)(minimum(v)) >= param.b0
    #     @warn "use higher order EOS!"
    # end
end

initparam(eos, ::Any) = fieldvalues(eos.param)
function initparam(eos::EnergyEoss, energies)
    if iszero(eos.param.e0)
        @set! eos.param.e0 = minimum(energies)
    end
    return fieldvalues(eos.param)
end

fieldvalues(x) = (getfield(x, i) for i in 1:nfields(x))  # Do not export!

end
