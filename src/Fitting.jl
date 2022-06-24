module Fitting

using ConstructionBase: constructorof, setproperties
using LsqFit: curve_fit, coef
using PolynomialRoots: roots
using Polynomials: fit, derivative, coeffs, derivative
using Unitful: AbstractQuantity, ustrip, unit, uconvert

using ..EquationsOfStateOfSolids:
    FiniteStrainParameters, Parameters, EnergyEquation, orderof
using ..FiniteStrains: FiniteStrain, ToStrain, FromStrain, Dⁿᵥf, straintype

export fiteos, linfit, nonlinfit

# See https://github.com/JuliaMath/Roots.jl/blob/bf0da62/src/utils.jl#L9-L11
struct ConvergenceFailed
    msg::String
end

struct CriterionNotMet
    msg::String
end

abstract type FittingMethod end
struct LinearFitting <: FittingMethod end
struct NonLinearFitting <: FittingMethod end

# ================================== Linear fitting ==================================
linfit(eos::EnergyEquation, volumes, energies; kwargs...) =
    EnergyEquation(fiteos(volumes, energies, eos.param, LinearFitting(); kwargs...))

"""
    fiteos(volumes, energies, initial_params::FiniteStrainParameters, LinearFitting(); kwargs...)

Fit an equation of state ``E(V)`` using linear algorithms.

# Arguments
- `maxiter::Integer=1000`: .
- `conv_thr::AbstractFloat=1e-12`: .
- `root_thr::AbstractFloat=1e-20`: .
- `verbose::Bool=false`: .

!!! note
    If you want to fit with `BigFloat` data, you need to install
    [`GenericSVD.jl`](https://github.com/JuliaLinearAlgebra/GenericSVD.jl) and `using GenericSVD`
    before fittting!
"""
function fiteos(
    volumes,
    energies,
    initial_params::FiniteStrainParameters,
    ::LinearFitting;
    maxiter = 1000,
    conv_thr = 1e-12,
    root_thr = 1e-20,
    verbose = false,
)
    deg, S = orderof(initial_params), straintype(initial_params)
    v0 = iszero(initial_params.v0) ? volumes[findmin(energies)[2]] : initial_params.v0  # Initial v0
    uv, ue = unit(v0), unit(energies[1])
    v0 = ustrip(v0)
    volumes = collect(map(x -> ustrip(uv, x), volumes))  # `parent` is needed to unwrap `DimArray`
    energies = collect(map(x -> ustrip(ue, x), energies))
    for i in 1:maxiter  # Self consistent loop
        strains = map(ToStrain{S}(v0), volumes)
        if !(isreal(strains) && isreal(energies))
            throw(DomainError("the strains or the energies are complex!"))
        end
        poly = fit(real(strains), real(energies), deg)
        f0, e0 = min_of_min(poly, root_thr)
        v0_prev, v0 = v0, FromStrain{S}(v0)(f0)  # Record v0 to v0_prev, then update v0
        if abs((v0_prev - v0) / v0_prev) <= conv_thr
            if verbose
                @info "convergence reached after $i steps!"
            end
            fᵥ = map(deg -> Dⁿᵥf(S(), deg, v0)(v0), 1:4)
            e_f = map(deg -> derivative(poly, deg)(f0), 1:4)
            b0, b′0, b″0 = Dₚb(fᵥ, e_f)
            param = update(
                initial_params;
                v0 = v0 * uv,
                b0 = b0(v0) * ue / uv,
                b′0 = b′0(v0),
                b″0 = b″0(v0) * uv / ue,
                e0 = e0 * ue,
            )
            return param
        end
    end
    throw(ConvergenceFailed("convergence not reached after $maxiter steps!"))
end

function islocalmin(x, y)  # `x` & `y` are both real
    y″ₓ = real(derivative(y, 2)(x))  # Must be real
    return y″ₓ > zero(y″ₓ)  # If 2nd derivative at x > 0, (x, y(x)) is a local minimum
end

function localmin(y, root_thr = 1e-20)  # `y` is a polynomial (could be complex)
    y′ = derivative(y, 1)
    pool = roots(coeffs(y′); polish = true, epsilon = root_thr)
    real_roots = real(filter(isreal, pool))  # Complex volumes are meaningless
    if isempty(real_roots)
        throw(CriterionNotMet("no real extrema found! Consider changing `root_thr`!"))  # For some polynomials, could be all complex
    else
        localminima = filter(x -> islocalmin(x, y), real_roots)
        if isempty(localminima)
            throw(CriterionNotMet("no local minima found!"))
        else
            return localminima
        end
    end
end

# https://stackoverflow.com/a/21367608/3260253
function min_of_min(y, root_thr = 1e-20)  # Find the minimum of the local minima
    localminima = localmin(y, root_thr)
    y0, i = findmin(map(y, localminima))  # `y0` must be real, or `findmap` will error
    x0 = localminima[i]
    return x0, y0
end

# See Eq. (55) - (57) in Ref. 1.
function Dₚb(fᵥ, e_f)  # Bulk modulus & its derivatives
    e″ᵥ = D²ᵥe(fᵥ, e_f)
    e‴ᵥ = D³ᵥe(fᵥ, e_f)
    b0 = v -> v * e″ᵥ
    b′0 = v -> -v * e‴ᵥ / e″ᵥ - 1
    b″0 = v -> (v * (D⁴ᵥe(fᵥ, e_f) * e″ᵥ - e‴ᵥ^2) + e‴ᵥ * e″ᵥ) / e″ᵥ^3
    return b0, b′0, b″0  # 3 lazy functions
end

function update(x::FiniteStrainParameters; kwargs...)
    patch = (; (f => kwargs[f] for f in propertynames(x))...)
    return setproperties(x, patch)
end

# Energy-volume derivatives, see Eq. (50) - (53) in Ref. 1.
D¹ᵥe(fᵥ, e_f) = e_f[1] * fᵥ[1]
D²ᵥe(fᵥ, e_f) = e_f[2] * fᵥ[1]^2 + e_f[1] * fᵥ[2]
D³ᵥe(fᵥ, e_f) = e_f[3] * fᵥ[1]^3 + 3fᵥ[1] * fᵥ[2] * e_f[2] + e_f[1] * fᵥ[3]
D⁴ᵥe(fᵥ, e_f) =
    e_f[4] * fᵥ[1]^4 +
    6fᵥ[1]^2 * fᵥ[2] * e_f[3] +
    (4fᵥ[1] * fᵥ[3] + 3fᵥ[3]^2) * e_f[2] +
    e_f[1] * fᵥ[4]

# ================================== Nonlinear fitting ==================================
nonlinfit(eos::EnergyEquation, volumes, energies; kwargs...) =
    EnergyEquation(fiteos(volumes, energies, eos.param, NonLinearFitting(); kwargs...))

"""
    nonlinfit(xs, ys, initial_params::Parameters, NonLinearFitting(); kwargs...)

Fit an equation of state ``E(V)`` using nonlinear algorithms.

# Arguments
- `xtol::AbstractFloat=1e-16`: .
- `gtol::AbstractFloat=1e-16`: .
- `maxiter::Integer=1000`: .
- `min_step_quality::AbstractFloat=1e-16`: .
- `good_step_quality::AbstractFloat=0.75`: .
- `verbose::Bool=false`: .
"""
function fiteos(
    volumes,
    energies,
    initial_params::T,
    ::NonLinearFitting;
    xtol = 1e-16,
    gtol = 1e-16,
    maxiter = 1000,
    min_step_quality = 1e-16,
    good_step_quality = 0.75,
    verbose = false,
) where {T<:Parameters}
    model = (volume, params) -> EnergyEquation(constructorof(T)(params...)).(volume)
    x, y, p0 = preprocess(volumes, energies, initial_params)
    fit = curve_fit(  # See https://github.com/JuliaNLSolvers/LsqFit.jl/blob/f687631/src/levenberg_marquardt.jl#L3-L28
        model,
        x,
        y,
        p0;
        x_tol = xtol,
        g_tol = gtol,
        maxIter = maxiter,
        min_step_quality = min_step_quality,
        good_step_quality = good_step_quality,
        show_trace = verbose,
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
        params = setproperties(params; e0 = uconvert(unit(params.e0), minimum(energies)))
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
            x |> up
        elseif f == :b″0
            x |> up^(-1)
        elseif f == :b‴0
            x |> up^(-2)
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
            x * up |> u
        elseif f == :b″0
            x * up^(-1) |> u
        elseif f == :b‴0
            x * up^(-2) |> u
        elseif f in (:v0, :b′0, :e0)
            x * u
        else
            error("unknown field `$f`!")
        end
    end
    return constructorof(typeof(p0))(params...)
end

end
