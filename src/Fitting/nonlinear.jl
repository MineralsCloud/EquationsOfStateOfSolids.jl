using ConstructionBase: constructorof
using LsqFit: curve_fit, coef

using ..EquationsOfStateOfSolids: Parameters

export nonlinfit

struct NonLinearFitting <: FittingMethod end

function nonlinfit(volumes, energies, initial_params; kwargs...)
    return fit(volumes, energies, initial_params, NonLinearFitting(); kwargs...)
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
    init_guess::Parameters,
    ::NonLinearFitting;
    xtol=1e-16,
    gtol=1e-16,
    maxiter=1000,
    min_step_quality=1e-16,
    good_step_quality=0.75,
    verbose=false,
)
    constructor = constructorof(typeof(init_guess))
    model = (volume, params) -> EnergyEquation(constructor(params...)).(volume)
    x, y, p0 = preprocess(volumes, energies, init_guess)
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
        params = reconstruct_units(coef(fit), init_guess)
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

unit(x) = oneunit(x) / one(x)

function preprocess(volumes, energies, params)  # Do not export!
    if iszero(params.e0)
        params = setproperties(params; e0=minimum(energies) / unit(params.e0))
    end
    energies ./= unit(params.e0)  # Unitless now
    volumes ./= unit(params.v0)  # Unitless now
    params = unitless(params)  # Unitless now
    return map(collect, (float.(volumes), float.(energies), float.(params)))
end

function unitless(params::Parameters)  # Normalize units of `params`
    punit = unit(params.e0) / unit(params.v0)  # Pressure/bulk modulus unit
    return Iterators.map(fieldnames(typeof(params))) do f
        x = getfield(params, f)
        if f == :b0
            x / punit
        elseif f == :b″0
            x * punit
        elseif f == :b‴0
            x * punit^2
        elseif f in (:v0, :b′0, :e0)
            x
        else
            error("unknown field `$f`!")
        end
    end
end

function reconstruct_units(p, init_guess)
    punit = unit(init_guess.e0) / unit(init_guess.v0)  # Pressure/bulk modulus unit
    params = Iterators.map(zip(p, fieldnames(typeof(init_guess)))) do (x, f)
        x0 = getfield(init_guess, f)
        unit(x0) * if f == :b0
            x * punit
        elseif f == :b″0
            x / punit
        elseif f == :b‴0e
            x / punit^2
        elseif f in (:v0, :b′0, :e0)
            x
        else
            error("unknown field `$f`!")
        end
    end
    return constructorof(typeof(init_guess))(params...)
end
