module Fitting

using ConstructionBase: constructorof, setproperties
using LsqFit: curve_fit, coef
using PolynomialRoots: roots
using Polynomials: Polynomial, fit, derivative, coeffs, derivative
using Serialization: serialize
using Unitful: AbstractQuantity, ustrip, unit, uconvert

using EquationsOfStateOfSolids: ispositive
using ..Collections:
    EquationOfStateOfSolids,
    FiniteStrainParameters,
    Parameters,
    PressureEOS,
    EnergyEOS,
    BulkModulusEOS,
    FiniteStrain,
    orderof,
    volume2strain,
    strain2volume,
    Dⁿᵥf,
    straintype

export linfit, nonlinfit, eosfit

# See https://github.com/JuliaMath/Roots.jl/blob/bf0da62/src/utils.jl#L9-L11
struct ConvergenceFailed
    msg::String
end

struct NoRootFound
    msg::String
end

eosfit(eos::EnergyEOS{<:FiniteStrainParameters}, volumes, energies; kwargs...) =
    linfit(eos, volumes, energies; kwargs...)
eosfit(eos, xs, ys; kwargs...) = nonlinfit(eos, xs, ys; kwargs...)

# ================================== Linear fitting ==================================
function linfit(
    eos::EnergyEOS{<:FiniteStrainParameters},
    volumes,
    energies;
    maxiter = 1000,
    conv_thr = 1e-12,
    root_thr = 1e-20,
    verbose = false,
)::FiniteStrainParameters
    deg = orderof(eos.param)
    s = straintype(eos.param)()
    v0 = iszero(eos.param.v0) ? volumes[findmin(energies)[2]] : eos.param.v0  # Initial v0
    uv, ue = unit(v0), unit(energies[1])
    uvrule = uconvert(unit(volumes[1]), 1 * uv)
    v0 = ustrip(v0)
    volumes = parent(map(x -> ustrip(uv, x), volumes))  # `parent` is needed to unwrap `DimArray`
    energies = parent(map(x -> ustrip(ue, x), energies))
    for i in 1:maxiter  # Self consistent loop
        strains = map(volume2strain(s, v0), volumes)
        poly = fit(strains, energies, deg)
        f0, e0 = _absminimum(poly, root_thr)
        v0_prev, v0 = v0, strain2volume(s, v0)(f0)  # Record v0 to v0_prev, then update v0
        if abs((v0_prev - v0) / v0_prev) <= conv_thr
            if verbose
                @info "convergence reached after $i steps!"
            end
            fᵥ = map(deg -> Dⁿᵥf(s, deg, v0)(v0), 1:4)
            e_f = map(deg -> derivative(poly, deg)(f0), 1:4)
            b0, b′0, b″0 = _Dₚb(fᵥ, e_f)
            return _update(
                eos.param;
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

function _islocalmin(x, y)
    y″ₓ = derivative(y, 2)(x)
    if isreal(y″ₓ)
        return ispositive(real(y″ₓ))  # If 2nd derivative at `x > 0`, `(x, y)` is a local minimum.
    else
        throw(DomainError("the 2nd derivative of the polynomial is a complex!"))
    end
end

function _localminima(y::Polynomial, root_thr = 1e-20)
    y′ = derivative(y, 1)
    rawpool = roots(coeffs(y′); polish = true, epsilon = root_thr)
    pool = real(filter(isreal, rawpool))  # Complex volumes are meaningless
    if isempty(pool)
        throw(NoRootFound("no real extrema found! Consider changing `root_thr`!"))  # For some polynomials, could be all complex
    else
        localminima = filter(x -> _islocalmin(x, y), pool)
        if isempty(localminima)
            throw(NoRootFound("no local minima found!"))
        else
            return localminima
        end
    end
end

# https://stackoverflow.com/a/21367608/3260253
function _absminimum(y, root_thr = 1e-20)  # Find the minimal in the minima
    localminima = _localminima(y, root_thr)
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
function nonlinfit(
    eos::EquationOfStateOfSolids,
    xs,
    ys;
    xtol = 1e-16,
    gtol = 1e-16,
    maxiter = 1000,
    min_step_quality = 1e-16,
    good_step_quality = 0.75,
    verbose = false,
)::Parameters
    model = createmodel(eos)
    p0, xs, ys = _prepare(eos, xs, ys)
    fit = curve_fit(  # See https://github.com/JuliaNLSolvers/LsqFit.jl/blob/f687631/src/levenberg_marquardt.jl#L3-L28
        model,
        xs,
        ys,
        collect(last.(p0));
        x_tol = xtol,
        g_tol = gtol,
        maxIter = maxiter,
        min_step_quality = min_step_quality,
        good_step_quality = good_step_quality,
        show_trace = verbose,
    )
    if fit.converged
        T = constructorof(typeof(eos.param))
        result = T((x * c for (x, c) in zip(coef(fit), first.(p0)))...)
        _checkresult(result)
        return result
    else
        throw(ConvergenceFailed("convergence not reached after $maxiter steps!"))
    end
end

function createmodel(::S) where {T,S<:EquationOfStateOfSolids{T}}  # Do not export!
    constructor = constructorof(S) ∘ constructorof(T)
    return (x, p) -> map(constructor(p), x)
end

function _checkresult(param::Parameters)  # Do not export!
    if param.v0 <= zero(param.v0) || param.b0 <= zero(param.b0)
        @error "either `v0 = $(param.v0)` or `b0 = $(param.b0)` is negative!"
    end
    # if PressureEoss(param)(minimum(v)) >= param.b0
    #     @warn "use higher order EOS!"
    # end
end

function _prepare(eos, xs, ys)  # Do not export!
    xs, ys = _collect_float(xs), _collect_float(ys)  # `xs` & `ys` may not be arrays
    if eos isa EnergyEOS && iszero(eos.param.e0)
        eos = EnergyEOS(setproperties(eos.param; e0 = minimum(ys)))  # Energy minimum as e0
    end
    return _ustrip_all(eos, xs, ys)
end

# No need to constrain `eltype`, `ustrip` will error if `Real` and `AbstractQuantity` are met.
function _ustrip_all(eos, xs, ys)  # Do not export!
    xs, ys = ustrip.(unit(eos.param.v0), xs), ustrip.(_yunit(eos), ys)
    punit = unit(eos.param.e0) / unit(eos.param.v0)
    return map(fieldnames(typeof(eos.param))) do f
        x = getfield(eos.param, f)
        u = unit(x)
        if f == :b0
            uconvert(u, 1 * punit) => ustrip(punit, float(x))
        elseif f == :b″0
            r = punit^(-1)
            uconvert(u, 1 * r) => ustrip(r, float(x))
        elseif f == :b‴0
            r = punit^(-2)
            uconvert(u, 1 * r) => ustrip(r, float(x))
        else
            oneunit(x) => ustrip(float(x))
        end
    end, xs, ys
end
_yunit(eos::EnergyEOS) = unit(eos.param.e0)
_yunit(eos::Union{PressureEOS,BulkModulusEOS}) = unit(eos.param.e0) / unit(eos.param.v0)

_collect_float(x) = collect(float.(x))  # Do not export!

_mapfields(f, x) = (f(getfield(x, i)) for i in 1:nfields(x))  # Do not export!

Base.float(p::Parameters) = constructorof(typeof(p))(_mapfields(float, p)...)  # Not used here but may be useful

end
