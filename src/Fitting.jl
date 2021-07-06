module Fitting

using ConstructionBase: constructorof, setproperties
using LsqFit: curve_fit, coef
using PolynomialRoots: roots
using Polynomials: fit, derivative, coeffs, derivative
using Unitful: AbstractQuantity, NoUnits, ustrip, unit, uconvert

using ..EquationsOfStateOfSolids:
    EquationOfStateOfSolids,
    FiniteStrainParameters,
    Parameters,
    PressureEquation,
    EnergyEquation,
    BulkModulusEquation,
    orderof,
    _ispositive,
    getparam
using ..FiniteStrains: FiniteStrain, ToStrain, FromStrain, Dⁿᵥf, straintype

import Unitful

export linfit, nonlinfit, eosfit

# See https://github.com/JuliaMath/Roots.jl/blob/bf0da62/src/utils.jl#L9-L11
struct ConvergenceFailed
    msg::String
end

struct CriterionNotMet
    msg::String
end

eosfit(eos::EnergyEquation{<:FiniteStrainParameters}, volumes, energies; kwargs...) =
    linfit(eos, volumes, energies; kwargs...)
eosfit(eos, xs, ys; kwargs...) = nonlinfit(eos, xs, ys; kwargs...)

# ================================== Linear fitting ==================================
"""
    linfit(eos::EnergyEquation{<:FiniteStrainParameters}, volumes, energies; kwargs...)

# Keyword Arguments
- `maxiter::Integer=1000`: .
- `conv_thr::AbstractFloat=1e-12`: .
- `root_thr::AbstractFloat=1e-20`: .
- `verbose::Bool=false`: .

!!! note
    If you want to fit with `BigFloat` data, you need to install
    [`GenericSVD.jl`](https://github.com/JuliaLinearAlgebra/GenericSVD.jl) and `using
    GenericSVD` before fittting!
"""
function linfit(
    eos::EnergyEquation{<:FiniteStrainParameters},
    volumes,
    energies;
    maxiter = 1000,
    conv_thr = 1e-12,
    root_thr = 1e-20,
    verbose = false,
)::FiniteStrainParameters
    deg = orderof(getparam(eos))
    S = straintype(getparam(eos))
    v0 = iszero(getparam(eos).v0) ? volumes[findmin(energies)[2]] : getparam(eos).v0  # Initial v0
    uv, ue = unit(v0), unit(energies[1])
    uvrule = uconvert(unit(volumes[1]), 1 * uv)
    v0 = ustrip(v0)
    volumes = collect(map(x -> ustrip(uv, x), volumes))  # `parent` is needed to unwrap `DimArray`
    energies = collect(map(x -> ustrip(ue, x), energies))
    for i in 1:maxiter  # Self consistent loop
        strains = map(ToStrain{S}(v0), volumes)
        if !(isreal(strains) && isreal(energies))
            throw(DomainError("the strains or the energies are complex!"))
        end
        poly = fit(real(strains), real(energies), deg)
        f0, e0 = _minofmin(poly, root_thr)
        v0_prev, v0 = v0, FromStrain{S}(v0)(f0)  # Record v0 to v0_prev, then update v0
        if abs((v0_prev - v0) / v0_prev) <= conv_thr
            if verbose
                @info "convergence reached after $i steps!"
            end
            fᵥ = map(deg -> Dⁿᵥf(S(), deg, v0)(v0), 1:4)
            e_f = map(deg -> derivative(poly, deg)(f0), 1:4)
            b0, b′0, b″0 = _Dₚb(fᵥ, e_f)
            return _update(
                getparam(eos);
                v0 = v0 * uv,
                b0 = b0(v0) * ue / uv,
                b′0 = b′0(v0),
                b″0 = b″0(v0) * uv / ue,
                e0 = e0 * ue,
            )
        end
    end
    throw(ConvergenceFailed("convergence not reached after $maxiter steps!"))
end

function _islocalmin(x, y)  # `x` & `y` are both real
    y″ₓ = derivative(y, 2)(x)  # Must be real
    return _ispositive(real(y″ₓ))  # If 2nd derivative at x > 0, (x, y(x)) is a local minimum
end

function _localmin(y, root_thr = 1e-20)  # `y` is a polynomial (could be complex)
    y′ = derivative(y, 1)
    pool = roots(coeffs(y′); polish = true, epsilon = root_thr)
    real_roots = real(filter(isreal, pool))  # Complex volumes are meaningless
    if isempty(real_roots)
        throw(CriterionNotMet("no real extrema found! Consider changing `root_thr`!"))  # For some polynomials, could be all complex
    else
        localminima = filter(x -> _islocalmin(x, y), real_roots)
        if isempty(localminima)
            throw(CriterionNotMet("no local minima found!"))
        else
            return localminima
        end
    end
end

# https://stackoverflow.com/a/21367608/3260253
function _minofmin(y, root_thr = 1e-20)  # Find the minimum of the local minima
    localminima = _localmin(y, root_thr)
    y0, i = findmin(map(y, localminima))  # `y0` must be real, or `findmap` will error
    x0 = localminima[i]
    return x0, y0
end

# See Eq. (55) - (57) in Ref. 1.
function _Dₚb(fᵥ, e_f)  # Bulk modulus & its derivatives
    e″ᵥ = _D²ᵥe(fᵥ, e_f)
    e‴ᵥ = _D³ᵥe(fᵥ, e_f)
    b0 = v -> v * e″ᵥ
    b′0 = v -> -v * e‴ᵥ / e″ᵥ - 1
    b″0 = v -> (v * (_D⁴ᵥe(fᵥ, e_f) * e″ᵥ - e‴ᵥ^2) + e‴ᵥ * e″ᵥ) / e″ᵥ^3
    return b0, b′0, b″0  # 3 lazy functions
end

# Energy-volume derivatives, see Eq. (50) - (53) in Ref. 1.
_D¹ᵥe(fᵥ, e_f) = e_f[1] * fᵥ[1]
_D²ᵥe(fᵥ, e_f) = e_f[2] * fᵥ[1]^2 + e_f[1] * fᵥ[2]
_D³ᵥe(fᵥ, e_f) = e_f[3] * fᵥ[1]^3 + 3fᵥ[1] * fᵥ[2] * e_f[2] + e_f[1] * fᵥ[3]
_D⁴ᵥe(fᵥ, e_f) =
    e_f[4] * fᵥ[1]^4 +
    6fᵥ[1]^2 * fᵥ[2] * e_f[3] +
    (4fᵥ[1] * fᵥ[3] + 3fᵥ[3]^2) * e_f[2] +
    e_f[1] * fᵥ[4]

function _update(x::FiniteStrainParameters; kwargs...)
    patch = (; (f => kwargs[f] for f in propertynames(x))...)
    return setproperties(x, patch)
end

# ================================== Nonlinear fitting ==================================
"""
    nonlinfit(eos::EquationOfStateOfSolids, xs, ys; kwargs...)

# Keyword Arguments
- `xtol::AbstractFloat=1e-16`: .
- `gtol::AbstractFloat=1e-16`: .
- `maxiter::Integer=1000`: .
- `min_step_quality::AbstractFloat=1e-16`: .
- `good_step_quality::AbstractFloat=0.75`: .
- `verbose::Bool=false`: .
"""
function nonlinfit(
    eos::EquationOfStateOfSolids,
    xdata,
    ydata;
    xtol = 1e-16,
    gtol = 1e-16,
    maxiter = 1000,
    min_step_quality = 1e-16,
    good_step_quality = 0.75,
    verbose = false,
)::Parameters
    model = buildmodel(eos)
    p0, xdata, ydata = _preprocess(eos, xdata, ydata)
    fit = curve_fit(  # See https://github.com/JuliaNLSolvers/LsqFit.jl/blob/f687631/src/levenberg_marquardt.jl#L3-L28
        model,
        xdata,
        ydata,
        p0;
        x_tol = xtol,
        g_tol = gtol,
        maxIter = maxiter,
        min_step_quality = min_step_quality,
        good_step_quality = good_step_quality,
        show_trace = verbose,
    )
    if fit.converged
        result = _postprocess(coef(fit), getparam(eos))
        checkresult(result)
        return result
    else
        throw(ConvergenceFailed("convergence not reached after $maxiter steps!"))
    end
end

# Do not export!
buildmodel(eos::EquationOfStateOfSolids{T}) where {T} =
    (x, p) -> constructorof(typeof(eos))(constructorof(T)(p...)).(x)

function checkresult(x::Parameters)  # Do not export!
    if x.v0 <= zero(x.v0) || x.b0 <= zero(x.b0)
        @error "either v0 ($(x.v0)) or b0 ($(x.b0)) is not positive!"
    end
    if PressureEquation(x)(x.v0) >= x.b0 && x isa FiniteStrainParameters
        @warn "consider using higher order EOS!"
    end
end

function _preprocess(eos, xdata, ydata)  # Do not export!
    p = getparam(eos)
    if eos isa EnergyEquation && iszero(p.e0)
        eos = EnergyEquation(setproperties(p; e0 = uconvert(unit(p.e0), minimum(ydata))))  # Energy minimum as e0, `uconvert` is important to keep the unit right!
    end
    return map(_float_collect, _unify(eos, xdata, ydata))  # `xs` & `ys` may not be arrays
end

function _unify(eos, xdata, ydata)  # Unify units of data
    p = getparam(eos)
    uy(eos::EnergyEquation) = unit(p.e0)
    uy(eos::Union{PressureEquation,BulkModulusEquation}) = unit(p.e0) / unit(p.v0)
    return ustrip.(_unormal(p)), ustrip.(unit(p.v0), xdata), ustrip.(uy(eos), ydata)
end

_unormal(p::Parameters) = collect(getfield(p, i) for i in 1:nfields(p))
function _unormal(p::Parameters{<:AbstractQuantity})  # Normalize units of `p`
    up = unit(p.e0) / unit(p.v0)  # Pressure/bulk modulus unit
    return map(fieldnames(typeof(p))) do f
        x = getfield(p, f)
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

_postprocess(p, p0::Parameters) = constructorof(typeof(p0))(p...)
function _postprocess(p, p0::Parameters{<:AbstractQuantity})
    up = unit(p0.e0) / unit(p0.v0)  # Pressure/bulk modulus unit
    param = map(enumerate(fieldnames(typeof(p0)))) do (i, f)
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
    return constructorof(typeof(p0))(param...)
end

_float_collect(x) = collect(float.(x))  # Do not export!

_fmap(f, x) = constructorof(typeof(x))((f(getfield(x, i)) for i in 1:nfields(x))...)  # Do not export!

"Convert all elements of a `Parameters` to floating point data types."
Base.float(p::Parameters) = _fmap(float, p)  # Not used here but may be useful

"Test whether all `p`'s elements are numerically equal to some real number."
Base.isreal(p::Parameters) = all(isreal(getfield(p, i)) for i in 1:nfields(p))  # Not used here but may be useful

"Construct a real `Parameters` from the real parts of the elements of p."
Base.real(p::Parameters) = _fmap(real, p)  # Not used here but may be useful

"""
    ustrip(p::Parameters)

Strip units from a `Parameters`.
"""
Unitful.ustrip(p::Parameters) = _fmap(ustrip, p)

end
