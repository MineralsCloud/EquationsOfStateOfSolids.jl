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
    p0, xs, ys = _preprocess(eos, xs, ys)
    fit = curve_fit(  # See https://github.com/JuliaNLSolvers/LsqFit.jl/blob/f687631/src/levenberg_marquardt.jl#L3-L28
        model,
        xs,
        ys,
        p0;
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

const collect_float = collect ∘ Base.Fix1(map, float)

function _preprocess(eos, xs, ys)  # Do not export!
    xs, ys = collect_float(xs), collect_float(ys)  # `xs` & `ys` may not be arrays
    eos, xs, ys = _unifyunit(eos, xs, ys)
    if eos isa EnergyEOS && iszero(eos.param.e0)
        @set! eos.param.e0 = minimum(ys)  # Energy minimum as e0
    end
    return collect(_splat(ustrip ∘ float, eos.param)), xs, ys
end

# Do not export!
function _unifyunit(eos::EnergyEOS{<:Parameters}, vs, es)
    vunit, eunit = unit(eos.param.v0), unit(eos.param.e0)
    punit = eunit / vunit
    vs, es = ustrip.(vunit, vs), ustrip.(eunit, es)
    if hasfield(typeof(eos.param), :b0)
        @set! eos.param.b0 = ustrip(punit, eos.param.b0)
    end
    if hasfield(typeof(eos.param), :b′′0)
        @set! eos.param.b0 = ustrip(punit^(-1), eos.param.b′′0)
    end
    if hasfield(typeof(eos.param), :b′′′0)
        @set! eos.param.b0 = ustrip(punit^(-2), eos.param.b′′′0)
    end
    return eos, vs, es
end
function _unifyunit(
    eos::Union{PressureEOS{T},BulkModulusEOS{T}},
    vs,
    ps,
) where {T<:Parameters}
    vunit, eunit = unit(eos.param.v0), unit(eos.param.e0)
    punit = eunit / vunit
    vs, ps = ustrip.(vunit, vs), ustrip.(punit, ps)
    if hasfield(typeof(eos.param), :b0)
        @set! eos.param.b0 = ustrip(punit, eos.param.b0)
    end
    if hasfield(typeof(eos.param), :b′′0)
        @set! eos.param.b0 = ustrip(punit^(-1), eos.param.b′′0)
    end
    if hasfield(typeof(eos.param), :b′′′0)
        @set! eos.param.b0 = ustrip(punit^(-2), eos.param.b′′′0)
    end
    return eos, vs, ps
end

_splat(x) = (getfield(x, i) for i in 1:nfields(x))  # Do not export!
_splat(f, x) = (f(getfield(x, i)) for i in 1:nfields(x))

Base.float(p::Parameters) = constructorof(typeof(p))(_splat(float, p)...)

end
