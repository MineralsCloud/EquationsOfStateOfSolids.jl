module Fitting

using ConstructionBase: constructorof
using LsqFit: curve_fit, coef
using Serialization: serialize
using Setfield: @set!
using Unitful: AbstractQuantity, ustrip, unit

using ..Collections:
    EquationOfStateOfSolids,
    FiniteStrainParameters,
    Parameters,
    PressureEOS,
    EnergyEOS,
    BulkModulusEOS

export linfit, nonlinfit

function nonlinfit(
    eos::EquationOfStateOfSolids{T},
    xs,
    ys;
    xtol = 1e-8,
    gtol = 1e-2,
    maxiter::Integer = 1000,
    min_step_quality = 1e-3,
    good_step_quality = 0.75,
    silent = true,
    saveto = "",
) where {T}
    model = createmodel(eos)
    p0 = initparam(eos, xs, ys)
    fit = curve_fit(  # See https://github.com/JuliaNLSolvers/LsqFit.jl/blob/f687631/src/levenberg_marquardt.jl#L3-L28
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
        _checkparam(param)
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

function _checkparam(param::FiniteStrainParameters)  # Do not export!
    if param.v0 <= zero(param.v0) || param.b0 <= zero(param.b0)
        @error "fitted `v0` or `b0` is negative!"
    end
    # if PressureEoss(param)(minimum(v)) >= param.b0
    #     @warn "use higher order EOS!"
    # end
end

function initparam(eos, xs, ys)
    _initparam(eos, xs, ys)
    return fieldvalues(eos.param)
end
function _initparam(eos::EnergyEOS, ::Any, energies)
    if iszero(eos.param.e0)
        @set! eos.param.e0 = minimum(energies)
    end
function _unifyunit(eos::EnergyEOS{<:Parameters{<:AbstractQuantity}}, vs, es)
    es = ustrip.(unit(eos.param.e0), es)
    vs = ustrip.(unit(eos.param.v0), vs)
    return vs, es
end
function _unifyunit(
    eos::Union{PressureEOS{T},BulkModulusEOS{T}},
    vs,
    ps,
) where {T<:Parameters{<:AbstractQuantity}}
    ps = ustrip.(unit(eos.param.e0) / unit(eos.param.v0), ps)
    vs = ustrip.(unit(eos.param.v0), vs)
    return vs, ps
end

_fieldvalues(x) = (getfield(x, i) for i in 1:nfields(x))  # Do not export!

end
