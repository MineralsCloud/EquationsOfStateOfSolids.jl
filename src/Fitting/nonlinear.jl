using ConstructionBase: constructorof
using LsqFit: curve_fit, coef
using Unitful: AbstractQuantity, uconvert

using ..EquationsOfStateOfSolids: Parameters

export nonlinfit

struct NonLinearFitting <: FittingMethod end

"""
    nonlinfit(eos::EnergyEquation, volumes, energies; kwargs...)

A wrapper for `fit` with nonlinear fitting method.

# Arguments
- `xtol::AbstractFloat=1e-16`: .
- `gtol::AbstractFloat=1e-16`: .
- `maxiter::Integer=1000`: .
- `min_step_quality::AbstractFloat=1e-16`: .
- `good_step_quality::AbstractFloat=0.75`: .
- `verbose::Bool=false`: .
"""
function nonlinfit(eos::EnergyEquation, volumes, energies; kwargs...)
    return EnergyEquation(fit(volumes, energies, eos.param, NonLinearFitting(); kwargs...))
end

"""
    fit(xs, ys, initial_params::Parameters, NonLinearFitting(); kwargs...)

Fit an equation of state ``E(V)`` using nonlinear algorithms.

# Arguments
- `xtol::AbstractFloat=1e-16`: .
- `gtol::AbstractFloat=1e-16`: .
- `maxiter::Integer=1000`: .
- `min_step_quality::AbstractFloat=1e-16`: .
- `good_step_quality::AbstractFloat=0.75`: .
- `verbose::Bool=false`: .
"""
function fit(
    volumes,
    energies,
    initial_params::T,
    ::NonLinearFitting;
    xtol=1e-16,
    gtol=1e-16,
    maxiter=1000,
    min_step_quality=1e-16,
    good_step_quality=0.75,
    verbose=false,
) where {T<:Parameters}
    model = (volume, params) -> EnergyEquation(constructorof(T)(params...)).(volume)
    x, y, p0 = preprocess(volumes, energies, initial_params)
    fit = curve_fit(  # See https://github.com/JuliaNLSolvers/LsqFit.jl/blob/f687631/src/levenberg_marquardt.jl#L3-L28
        model,
        x,
        y,
        p0;
        x_tol=xtol,
        g_tol=gtol,
        maxIter=maxiter,
        min_step_quality=min_step_quality,
        good_step_quality=good_step_quality,
        show_trace=verbose,
    )
    if fit.converged
        params = reconstruct_params(coef(fit), initial_params)
        checkresult(params)
        return params
    else
        throw(ConvergenceFailed("convergence not reached after $maxiter steps!"))
    end
end

function checkresult(params::Parameters)  # Do not export!
    if params.v0 <= zero(params.v0)
        @info "V₀ ($(params.v0)) is not positive!"
    end
    if params.b0 <= zero(params.b0)
        @info "B₀ ($(params.b0)) is not positive!"
    end
end

function preprocess(volumes, energies, params)  # Do not export!
    volumes = ustrip.(unit(params.v0), volumes)  # Unify units of data
    if iszero(params.e0)
        # Energy minimum as e0, `uconvert` is important to keep the unit right!
        params = setproperties(params; e0=uconvert(unit(params.e0), minimum(energies)))
    end
    energies = ustrip.(unit(params.e0), energies)
    params = ustrip.(unormalize(params))
    return map(collect, (float.(volumes), float.(energies), float.(params)))
end

unormalize(params::Parameters) = (getfield(params, f) for f in fieldnames(typeof(params)))
function unormalize(params::Parameters{<:AbstractQuantity})  # Normalize units of `params`
    up = unit(params.e0) / unit(params.v0)  # Pressure/bulk modulus unit
    return Iterators.map(fieldnames(typeof(params))) do f
        x = getfield(params, f)
        if f == :b0
            up(x)
        elseif f == :b″0
            up^(-1)(x)
        elseif f == :b‴0
            up^(-2)(x)
        elseif f in (:v0, :b′0, :e0)
            x
        else
            error("unknown field `$f`!")
        end
    end
end

reconstruct_params(p, p0::Parameters) = constructorof(typeof(p0))(p...)
function reconstruct_params(p, p0::Parameters{<:AbstractQuantity})
    up = unit(p0.e0) / unit(p0.v0)  # Pressure/bulk modulus unit
    params = Iterators.map(enumerate(fieldnames(typeof(p0)))) do (i, f)
        x = p[i]
        u = unit(getfield(p0, f))
        if f == :b0
            u(x * up)
        elseif f == :b″0
            u(x * up^(-1))
        elseif f == :b‴0
            u(x * up^(-2))
        elseif f in (:v0, :b′0, :e0)
            x * u
        else
            error("unknown field `$f`!")
        end
    end
    return constructorof(typeof(p0))(params...)
end
